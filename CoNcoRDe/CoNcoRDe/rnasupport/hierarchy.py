import random
import pickle
from copy import deepcopy
import RNA

COMPLEMENT = {
  "G": "UC",
  "A": "U",
  "U": "GA",
  "C": "G"
}

def randcomp(seq):
  return "".join([
    random.choice(COMPLEMENT[item])
    for item in reversed(seq)
  ])

class Block:
  def __init__(self, kind, size, children=None):
    self.kind = kind
    self.size = size
    self.capacity = 0
    self.sequence = None
    self.structure = None
    self.matches = False
    self.children = children or []

  def say(self):
    if self.kind == "S":
      children = "".join(map(lambda x: x.say(), self.children))
      if not children:
        children = "&"
      return "(" * self.size + children + ")" * self.size
    if self.kind == "L":
      return "." * self.size
    if self.kind in ("O", "C"):
      return "".join(map(lambda x: x.say(), self.children))

  def is_terminal(self):
    result = self.kind == "S" and len(self.children) == 1 and self.children[0].kind == "L"
    result = result or self.kind == "L"
    return result

  def as_tuple(self):
    return (self.kind, self.size, tuple(map(lambda x: x.as_tuple(), self.children)))

  def append(self, child):
    self.children.append(child)

  def set_sequence(self, seq, start=0):
    if self.kind in ("O", "C"):
      for child in self.children:
        start = child.set_sequence(seq, start=start)
      return start
    if self.kind == "L":
      stop = start + self.size
      self.sequence = seq[start:stop]
      start = stop
      return start
    if self.kind == "S":
      first_stop = start + self.size
      first_sequence = seq[start:first_stop]
      start = first_stop
      for child in self.children:
        start = child.set_sequence(seq, start=start)
      second_stop = start + self.size
      second_sequence = seq[start:second_stop]
      self.sequence = (first_sequence, second_sequence)
      start = second_stop
      return start

  def set_structure(self, seq, start=0):
    if self.kind in ("O", "C"):
      for child in self.children:
        start = child.set_structure(seq, start=start)
      return start
    if self.kind == "L":
      stop = start + self.size
      self.structure = seq[start:stop]
      start = stop
      return start
    if self.kind == "S":
      first_stop = start + self.size
      first_sequence = seq[start:first_stop]
      start = first_stop
      for child in self.children:
        start = child.set_structure(seq, start=start)
      second_stop = start + self.size
      second_sequence = seq[start:second_stop]
      self.structure = (first_sequence, second_sequence)
      start = second_stop
      return start

  def get_structure(self, recursive=True):
    if self.kind == "L":
      return self.structure
    if self.kind == "S":
      subseq = "&"
      if recursive:
        subseq = "".join([
          child.get_structure(recursive=True)
          for child in self.children
        ])
      return ("{}" + subseq + "{}").format(*self.structure)
    if self.kind in ("O", "C"):
      return "".join([
        child.get_structure(recursive=True)
        for child in self.children
      ])

  def match_structure(self):
    self.matches = self.say() == self.get_structure()
    for child in self.children:
      child.match_structure()

  def set_random_sequence(self, recursive=True):
    if recursive and self.kind in ("O", "C"):
      for child in self.children:
        child = child.set_random_sequence(recursive=True)
    if self.kind == "L":
      self.sequence = "".join(random.choices("GUAC", k=self.size))
    if self.kind == "S":
      sequence = "".join(random.choices("GUAC", k=self.size))
      self.sequence = (
        sequence,
        randcomp(sequence)
      )
      if recursive:
        for child in self.children:
          child = child.set_random_sequence(recursive=True)
    return self

  def mutate_sequence(self, recursive=True):
    if recursive and self.kind in ("O", "C"):
      for child in self.children:
        child = child.set_random_sequence(recursive=True)
    if self.kind == "L":
      self.sequence = "".join([
        random.choice("GUAC") if random.random() < 0.1 else item
        for item in self.sequence
      ])
    if self.kind == "S":
      sequence = "".join([
        random.choice("GUAC") if random.random() < 0.1 else item
        for item in self.sequence[0]
      ])
      self.sequence = (
        sequence,
        randcomp(sequence)
      )
      if recursive:
        for child in self.children:
          child = child.set_random_sequence(recursive=True)
    return self

  def get_sequence(self, recursive=True):
    if self.kind == "L":
      return self.sequence
    if self.kind == "S":
      subseq = "&"
      if recursive:
        subseq = "".join([
          child.get_sequence(recursive=True)
          for child in self.children
        ])
      return ("{}" + subseq + "{}").format(*self.sequence)
    if self.kind in ("O", "C"):
      return "".join([
        child.get_sequence(recursive=True)
        for child in self.children
      ])

  def cut(self, depth):
    result = Block(self.kind, self.size, [])
    result.sequence = self.sequence
    if depth > 0:
      for child in self.children:
        result.append(child.cut(depth - 1))
    return result

  def partial_substructures(self, depth=2):
    for substructure in self.substructures():
      yield substructure.cut(depth)

  def substructures(self):
    if not self.kind in ("O", "C"):
      yield self
    for child in self.children:
      if child.kind != "L":
        for substructure in child.substructures():
          yield substructure

  def mix(self, other, pswitch=0.1):
    new_children = [
      deepcopy(c2.mix(c1, pswitch=1 - pswitch)) if random.random() < pswitch else c1
      for c1, c2 in zip(self.children, other.children)
    ]
    self.children = new_children
    return self

  def random_substructure(self):
    substructures = list(self.substructures())
    return random.choice(substructures[1:])

  def random_component(self):
    if self.kind == "C":
      return random.choice(self.children)
    if self.kind == "L":
      raise ValueError("Has no components!")
    if len(self.children) == 1:
      return self.children[0]
    left = random.randrange(0, len(self.children))
    right = random.randrange(left + 1, len(self.children) + 1)
    subset = self.children[left:right]
    result = component(subset)
    if len(subset) == 1:
      result = subset[0]
    result = random.choice(subset)
    return result

  def __len__(self):
    if self.kind == "L":
      return self.size
    if self.kind == "S":
      return 2 * self.size + sum(map(len, self.children))
    if self.kind in ("O", "C"):
      return sum(map(len, self.children))

def stem(size, children=None):
  if children:
    return Block("S", size, children)
  return Block("S", size)

def loop(size):
  return Block("L", size)

def origin():
  return Block("O", 0)

def component(children):
  return Block("C", 0, children)

def _dot_mode(point, current_block, stack):
  mode = "."
  if point == ".":
    current_block.size += 1
  if point == "(":
    new_current = stack.pop()
    new_current.append(current_block)
    current_block = new_current
    stack.append(current_block)
    current_block = stem(0)
    current_block.capacity += 1
    mode = "("
  if point == ")":
    new_current = stack.pop()
    new_current.append(current_block)
    current_block = new_current
    current_block.size += 1
    current_block.capacity -= 1
    mode = ")"
  return mode, current_block

def _open_mode(point, current_block, stack):
  mode = "("
  if point == ".":
    mode = "."
    stack.append(current_block)
    current_block = loop(1)
  if point == "(":
    current_block.capacity += 1
  return mode, current_block

def _close_mode(point, current_block, stack):
  mode = ")"
  if point == ".":
    if current_block.capacity != 0:
      cap = current_block.capacity
      current_block = stem(0, children=[current_block])
      current_block.capacity = cap
      stack.append(current_block)
      current_block = loop(1)
    else:
      new_current = stack.pop()
      new_current.append(current_block)
      stack.append(new_current)
      current_block = loop(1)
    mode = "."
  if point == ")":
    if current_block.capacity == 0:
      new_current = stack.pop()
      new_current.append(current_block)
      current_block = new_current
    current_block.capacity -= 1
    current_block.size += 1
  if point == "(":
    if current_block.capacity != 0:
      cap = current_block.capacity
      current_block = stem(0, children=[current_block])
      current_block.capacity = cap
      stack.append(current_block)
      current_block = stem(0)
      current_block.capacity += 1
    else:
      new_current = stack.pop()
      new_current.append(current_block)
      stack.append(new_current)
      current_block = stem(0)
      current_block.capacity += 1
    mode = "("
  return mode, current_block

def parse_hierarchy(struct):
  stack = []
  current_block = origin()
  mode = "#"
  for point in struct:
    if mode == "#":
      stack.append(current_block)
      if point == ".":
        current_block = loop(0)
        mode = "."
      if point == "(":
        current_block = stem(0)
        mode = "("
    if mode == ".":
      mode, current_block = _dot_mode(point, current_block, stack)
    elif mode == "(":
      mode, current_block = _open_mode(point, current_block, stack)
    elif mode == ")":
      mode, current_block = _close_mode(point, current_block, stack)
  if stack:
    new_current = stack.pop()
    new_current.append(current_block)
    current_block = new_current
  return current_block

class HierarchyDesign:
  def __init__(self, population_size=10, library_regen=0.5,
               steps=100, inner_steps=5, inner_population_size=5):
    self.library = {}
    self.library_regen = library_regen
    self.population_size = population_size
    self.steps = steps
    self.inner_steps = inner_steps
    self.inner_population_size = inner_population_size

  def dump_library(self, path):
    with open(path, "wb") as f:
      pickle.dump(self.library, f)

  def load_library(self, path):
    with open(path, "rb") as f:
      self.library = pickle.load(f)

  def randomise(self, structure):
    structure.set_random_sequence()
    return structure

  def store(self, structure, sequence):
    if structure not in self.library:
      self.library[structure] = []
    if sequence not in self.library[structure]:
      self.library[structure].append(sequence)

  def init_library(self, structure):
    structure.set_random_sequence()
    for sub in structure.substructures():
      if sub.kind != "L" and self.is_designable(sub) and not sub.is_terminal():
        seq = sub.get_sequence()
        dot = sub.say()
        if len(seq) < 32:
          for idx in range(100):
            seq, error = RNA.inverse_fold(
              "".join(random.choices("GAUC", k=len(seq))),
              dot
            )
            struc, _ = RNA.fold(seq)
            if error == 0.0:
              print("solved substructure", dot, seq)
              self.store(struc, seq)

  def is_designable(self, structure):
    if structure.kind == "S" and structure.size < 2:
      return False
    if structure.kind == "C":
      for child in structure.children:
        if not self.is_designable(child):
          return False
    return True

  def expand_library(self, structure):
    substructures = []
    for sub in structure.substructures():
      dot = sub.say()
      count = len(self.library[dot]) if dot in self.library else 0
      substructures.append((dot, sub, count))
    max_count = max(map(lambda x: x[2], substructures))
    substructures = [
      x
      for sub in substructures
      for x in [sub[1]] * (max_count - sub[2] + 1)
    ]
    sub = random.choice(substructures)
    if not self.is_designable(sub):
      sub = deepcopy(sub)
      sub.mutate_sequence()
      sequence = sub.get_sequence()
      structure, _ = RNA.fold(sequence)
    else:
      sequence, structure = self.design(
        sub,
        steps=self.inner_steps,
        population_size=self.inner_population_size
      )
    self.store(structure, sequence)

  def randomise_substructure(self, structure):
    if structure.kind == "S":
      if random.random() < 0.1:
        structure.set_random_sequence(recursive=False)
    try:
      sub = structure.random_component()
    except Exception as e:
      return structure
    if sub.kind == "L":
      sub.set_random_sequence()
      return structure
    if random.random() < 0.1:
      sub.mutate_sequence()
      return structure

    self.from_library(sub)
    return structure

  def from_library(self, sub, design=True):
    dot = sub.say()
    if not dot in self.library:
      if design and self.is_designable(sub):
        seq = sub.get_sequence()
        if len(seq) < 64:
          seq, error = RNA.inverse_fold(seq, dot)
          struc, _ = RNA.fold(seq)
          self.store(struc, seq)
          if error == 0.0:
            print("solved substructure", struc, seq)
            sub.set_sequence(seq)
            return
        sequence, structure = self.design(
          sub,
          steps=self.inner_steps,
          population_size=self.inner_population_size
        )
        if structure == dot:
          print("solved substructure", structure, sequence)
          sub.set_sequence(sequence)
        else:
          sub.set_random_sequence()
      else:
        sub.set_random_sequence()
    else:
      sequence = random.choice(self.library[dot])
      sub.set_sequence(sequence)

    if random.random() < self.library_regen:
      if design and self.is_designable(sub):
        sequence, structure = self.design(
          sub,
          steps=self.inner_steps,
          population_size=self.inner_population_size
        )
      else:
        sub = deepcopy(sub)
        sub.mutate_sequence()
        sequence = sub.get_sequence()
        structure, _ = RNA.fold(sequence)
      dot = sub.say()
      if dot not in self.library and dot == structure:
        print("solved substructure", structure, sequence)
      self.store(structure, sequence)

  def fitness(self, structure):
    struc = structure.say()
    seq = structure.get_sequence()
    real_struc, E = RNA.fold(seq)
    mismatches = sum(int(x != y) for x, y in zip(struc, real_struc))
    return 100 * mismatches + E, seq, real_struc

  def assess(self, structure):
    fit, seq, struc = self.fitness(structure)
    self.store(struc, seq)
    return fit, seq, struc, structure

  def mutate(self, structure):
    result = deepcopy(structure)
    for idx in range(random.randrange(0, 10)):
      self.randomise_substructure(result)
    return result

  def design(self, structure, steps=None, population_size=None,
             early=False):
    upper = steps is None
    steps = steps or self.steps
    population_size = population_size or self.population_size
    population = [
      deepcopy(structure)
      for idx in range(population_size)
    ]
    for structure in population:
      self.from_library(structure, design=False)
    best = None
    best_struc = None
    best_fit = None
    for idx in range(steps):
      if upper:
        print("scoring...")
      fitness = sorted([
        self.assess(structure)
        for structure in population
      ], key=lambda x: x[0])
      fit, seq, struc, _ = fitness[0]
      if not best or fit < best_fit:
        best_fit = fit
        best = seq
        best_struc = struc
      if upper and fit < 0:
        print(seq, struc)
        if early:
          return seq, struc
      breed = [
        fitness[idx][3]
        for idx in range(5)
      ]
      if upper:
        print("mutating...")
      population = [
        self.mutate(random.choice(breed))
        for idx in range(len(population))
      ]
    return best, best_struc
