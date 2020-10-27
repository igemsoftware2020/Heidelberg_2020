import random
from RNA import fold, eval_structure_simple

import torch

from torchsupport.modules import one_hot_encode
from torchsupport.interacting.environments.environment import Environment
from torchsupport.data.namedtuple import namedtuple

from rnasupport.secstruct import block_bonds, parse_secondary_structure

class BlockRNA(Environment):
  def __init__(self, min_size=10, max_size=200, max_steps=1000):
    self.min_size = min_size
    self.max_size = max_size
    self.initialized = False
    self.sequence = None
    self.target = None
    self.state = None
    self.done = False
    self.max_steps = max_steps
    self.steps = 0

  def reset(self):
    size = random.randrange(self.min_size, self.max_size)
    seq = random.choices(["G", "A", "U", "C"], k=size)
    sequence = "".join(seq)
    structure, _ = fold(sequence)
    start = sequence
    start_structure = structure
    seq = random.choices(["G", "A", "U", "C"], k=size)
    sequence = "".join(seq)
    structure, _ = fold(sequence)
    target = structure

    mask = torch.zeros(self.max_size, dtype=torch.bool)
    mask[:size] = True
    target_bonds = torch.zeros(self.max_size, dtype=torch.long)
    b, _ = block_bonds(target)
    target_bonds[:size] = b
    start_bonds = torch.zeros(self.max_size, dtype=torch.long)
    b, _ = block_bonds(start_structure)
    start_bonds[:size] = b
    sequence_onehot = torch.zeros(4, self.max_size)
    sequence_onehot[:, :size] = one_hot_encode(start, "GAUC")

    self.target = target
    self.state = (
      sequence_onehot, start_bonds, target_bonds, mask
    )

    self.initialized = True
    self.done = False
    self.steps = 0

  def act(self, action):
    # update base using action:
    self.state[0][:, action[1]] = 0
    self.state[0][action[0], action[1]] = 1
    new_sequence = "".join([
      "GAUC"[item]
      for item in self.state[0].argmax(dim=0)
    ])

    # fold new structure
    structure, _ = fold(new_sequence)
    new_bonds = torch.zeros(self.max_size, dtype=torch.long)
    b, _ = block_bonds(structure)
    new_bonds[:len(self.sequence)] = b

    # compute new reward
    reward = -sum([
      int(x != y)
      for x, y in zip(structure, self.target)
    ])
    self.state[1] = new_bonds

    self.done = self.steps > self.max_steps
    return torch.tensor([reward])

  def observe(self):
    return self.state

  def is_done(self):
    return self.done

  @property
  def action_space(self):
    return None

  @property
  def observation_space(self):
    return None

  def schema(self):
    state = (
      torch.zeros(4, self.max_size),
      torch.zeros(self.max_size, dtype=torch.long),
      torch.zeros(self.max_size, dtype=torch.long),
      torch.zeros(self.max_size, dtype=torch.bool)
    )
    reward = torch.tensor(0.0)
    done = torch.tensor(0)
    action = torch.tensor([0, 0])
    sample = self.data_type(
      state=state, action=action, rewards=reward, done=done
    )
    return sample

class SeqRNA(BlockRNA):
  data_type = namedtuple("Data", ["state", "action", "rewards", "done"])
  state_type = namedtuple("State", ["sequence", "structure", "target", "mask", "position"])
  def __init__(self, min_size=10, max_size=200, max_steps=1000, close=False,
               weight_reward=True, norm_reward=True):
    self.min_size = min_size
    self.max_size = max_size
    self.close = close
    self.initialized = False
    self.sequence = None
    self.target = None
    self.structure = None
    self.pairs = None
    self.done = False
    self.count = 0
    self.place = 0

    self.memory = {}
    self.unsolved = []
    self.max_seen = 0
    self.weight_reward = weight_reward
    self.norm_reward = norm_reward

  def store(self, structure, solved=False):
    if structure not in self.memory:
      self.memory[structure] = [0, 0]
    self.memory[structure][int(solved)] += 1
    if self.memory[structure][0] > self.max_seen:
      self.max_seen = self.memory[structure][0]

  def solve_reward(self, structure):
    if self.weight_reward:
      seen, solved = self.memory[structure]
      return 100 * seen / (solved + 1)
    else:
      result = 100 if self.norm_reward else len(structure)
      return result

  def pack(self, sequence, current, target, position):
    size = len(sequence)
    mask = torch.zeros(self.max_size, dtype=torch.bool)
    mask[:size] = True
    target_bonds = torch.zeros(self.max_size, dtype=torch.long)
    b, _ = block_bonds(target)
    target_bonds[:size] = b
    start_bonds = torch.zeros(self.max_size, dtype=torch.long)
    b, _ = block_bonds(current)
    start_bonds[:size] = b
    sequence_onehot = torch.zeros(5, self.max_size)
    sequence_onehot[:, :size] = one_hot_encode(sequence, "GAUCN")
    position = torch.tensor([position], dtype=torch.long)
    return self.state_type(
      sequence=sequence_onehot,
      structure=start_bonds,
      target=target_bonds,
      mask=mask,
      position=position
    )

  def reset(self):
    size = random.randrange(self.min_size, self.max_size)
    seq = random.choices(["G", "A", "U", "C"], k=size)
    sequence = "".join(seq)
    structure, _ = fold(sequence)
    target = structure

    if self.unsolved and random.random() < 0.2:
      random.shuffle(self.unsolved)
      target = self.unsolved.pop()
    
    size = len(target)
    seq = "N" * size
    sequence = "".join(seq)
    structure, _ = fold(sequence)
    start = sequence
    start_structure = structure
    
    self.target = target
    self.store(target)
    self.sequence = start
    self.structure = start_structure
    self.pairs, *_ = parse_secondary_structure(self.target)
    self.place = 0
    self.count = 0

    self.initialized = True
    self.done = False
    self.steps = 0

  def hamming(self, struc, target):
    result = 0.0
    for s, t in zip(struc, target):
      if s != t:
        result += 1
    return 1 - result / len(struc)

  def complement(self, base):
    comp = dict(A="U", U="A", G="C", C="G")
    return comp[base]

  def act(self, action):
    # update base using action:
    reward = 0.0
    self.sequence = self.sequence[:self.place] + "GAUC"[action] + self.sequence[self.place + 1:]
    if self.close and self.place in self.pairs:
      pos = self.pairs[self.place]
      comp = self.complement(self.sequence[self.place])
      self.sequence = self.sequence[:pos] + comp + self.sequence[pos + 1:]
    self.structure, _ = fold(self.sequence)

    self.place = self.sequence.find("N")
    if self.place == -1:
      hamming = self.hamming(self.structure, self.target)
      reward += hamming
      self.count += 1
      if self.structure == self.target:
        reward = self.solve_reward(self.target)
        self.store(self.target, solved=True)
        self.done = True
      if self.count == 10:
        self.unsolved.append(self.target)
        self.done = True
      if self.place == -1:
        self.sequence = "".join([
          "N" if struc != target else pos
          for pos, struc, target in zip(
            self.sequence, self.structure, self.target
          )
        ])
        self.place = self.sequence.find("N")
      self.structure, _ = fold(self.sequence)

    return torch.tensor([reward])

  def observe(self):
    return self.pack(
      self.sequence, self.structure,
      self.target, self.place
    )

  def is_done(self):
    return self.done

  @property
  def action_space(self):
    return None

  @property
  def observation_space(self):
    return None

  def schema(self):
    state = self.state_type(
      sequence=torch.zeros(5, self.max_size),
      structure=torch.zeros(self.max_size, dtype=torch.long),
      target=torch.zeros(self.max_size, dtype=torch.long),
      mask=torch.zeros(self.max_size, dtype=torch.bool),
      position=torch.zeros(1, dtype=torch.long)
    )
    reward = torch.tensor(0.0)
    done = torch.tensor(0)
    action = torch.tensor(0)
    sample = self.data_type(
      state=state, action=action, rewards=reward, done=done
    )
    return sample

class PosRNA(SeqRNA):
  def reset(self):
    size = random.randrange(self.min_size, self.max_size)
    seq = random.choices(["G", "A", "U", "C"], k=size)
    sequence = "".join(seq)
    structure, _ = fold(sequence)
    target = structure

    if self.unsolved and random.random() < 0.2:
      random.shuffle(self.unsolved)
      target = self.unsolved.pop()
    size = len(target)

    self.target = target
    self.store(self.target)
    self.pairs, *_ = parse_secondary_structure(self.target)

    seq = "N" * size
    start = ""
    for idx in range(size):
      if idx in self.pairs and self.pairs[idx] < idx:
        start += self.complement(start[self.pairs[idx]])
      else:
        start += random.choice("GAUC")
    start = "".join(start)
    start_structure, _ = fold(start)

    self.sequence = start
    self.structure = start_structure
    self.place = 0
    self.count = 0

    self.previous = 0

    self.initialized = True
    self.done = False
    self.steps = 0

  def act(self, action):
    closing = {"A": "U", "U": "GA", "C": "G", "G": "UC"}
    # update base using action:
    reward = 0.0
    place = action % self.max_size
    choice = action // self.max_size
    self.sequence = self.sequence[:place] + "GAUC"[choice] + self.sequence[place + 1:]
    if self.close and place in self.pairs:
      pos = self.pairs[place]
      if self.sequence[pos] not in closing[self.sequence[place]]:
        comp = self.complement(self.sequence[pos])
        self.sequence = self.sequence[:pos] + comp + self.sequence[pos + 1:]
    self.structure, _ = fold(self.sequence)

    hamming = self.hamming(self.structure, self.target)
    reward = hamming - self.previous
    self.previous = hamming
    self.count += 1
    if self.structure == self.target:
      reward = self.solve_reward(self.target)
      self.store(self.target, solved=True)
      self.done = True
    if self.count == 2 * self.max_size:
      self.unsolved.append(self.target)
      self.done = True
    self.structure, _ = fold(self.sequence)

    return torch.tensor([reward])
