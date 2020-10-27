import os
import random

from RNA import fold

import torch
from torch.utils.data import Dataset

from torchsupport.modules.basic import one_hot_encode
from torchsupport.structured import PackedTensor, SubgraphStructure

from rnasupport.secstruct import tensor_structure, block_bonds

class RefoldData(Dataset):
  def __init__(self, path, mode="Train"):
    self.path = os.path.join(path, mode)
    self.structure_map = self.read_structure_map(self.path)
    self.sizes = [key for key in self.structure_map]

  def read_structure_map(self, path):
    result = {}
    for file_path in os.listdir(path):
      full_path = os.path.join(path, file_path)
      if os.path.isfile(full_path) and full_path.endswith("csv"):
        size = int(file_path.split(".")[0].split("-")[1])
        if size not in result:
          result[size] = []
        with open(full_path) as f:
          for line in f:
            if line.strip():
              source, source_structure, E_s, target, target_structure, E_t = line.strip().split(",")
              E_s = float(E_s)
              E_t = float(E_t)
              result[size].append((source, source_structure, E_s, target, target_structure, E_t))
      return result

  def __getitem__(self, index):
    size = random.choice(self.sizes)
    pick = random.choice(self.structure_map[size])
    return pick

  def __len__(self):
    return 10000

class FoldData(RefoldData):
  def read_structure_map(self, path):
    old_result = super().read_structure_map(path)
    result = {
      key : [
        item
        for seq_0, struc_0, E_0, seq_1, struc_1, E_1 in old_result[key]
        for item in [(seq_0, struc_0, E_0), (seq_1, struc_1, E_1)]
      ]
      for key in old_result
    }
    return result

class RNAGraphFoldData(FoldData):
  def __getitem__(self, index):
    sequence, structure, _ = super().__getitem__(index)

    sequence_gt = one_hot_encode(sequence, "GAUC", numeric=True)
    sequence_onehot = one_hot_encode(sequence, "GAUC").permute(1, 0)
    mask = torch.rand(len(sequence)) < random.random()
    sequence_onehot[mask] = 0
    sequence_onehot = torch.cat((sequence_onehot, mask[:, None].float()), dim=1)
    bond, stem, unbound = tensor_structure(structure)
    rna = SubgraphStructure(torch.zeros(len(sequence), dtype=torch.long))

    inputs = (
      PackedTensor(sequence_onehot), rna, bond, stem, unbound
    )
    outputs = PackedTensor(sequence_gt, split=False)

    return inputs, (outputs, PackedTensor(mask, split=False))

class RNAGraphFoldRecomputeData(Dataset):
  def __getitem__(self, index):
    size = random.randrange(10, 200)
    seq = random.choices(["G", "A", "U", "C"], k=size)
    sequence = "".join(seq)
    structure, _ = fold(sequence)

    sequence_gt = one_hot_encode(sequence, "GAUC", numeric=True)
    sequence_onehot = one_hot_encode(sequence, "GAUC").permute(1, 0)
    mask = torch.rand(len(sequence)) < random.random()
    if mask.sum() == 0:
      mask[random.randrange(len(sequence))] = 1
    sequence_onehot[mask] = 0
    sequence_onehot = torch.cat((sequence_onehot, mask[:, None].float()), dim=1)
    bond, stem, unbound = tensor_structure(structure)
    rna = SubgraphStructure(torch.zeros(len(sequence), dtype=torch.long))

    inputs = (
      PackedTensor(sequence_onehot), rna, bond, stem, unbound
    )
    outputs = PackedTensor(sequence_gt, split=False)

    return inputs, (outputs, PackedTensor(mask, split=False))

  def __len__(self):
    return 100000

class RNABlockFoldRecomputeData(RNAGraphFoldRecomputeData):
  def __init__(self, min_size=10, max_size=200):
    self.min_size = min_size
    self.max_size = max_size

  def __getitem__(self, index):
    size = random.randrange(self.min_size, self.max_size)
    seq = random.choices(["G", "A", "U", "C"], k=size)
    sequence = "".join(seq)
    structure, _ = fold(sequence)

    sequence_gt = torch.zeros(self.max_size, dtype=torch.long)
    sequence_gt[:size] = one_hot_encode(sequence, "GAUC", numeric=True)
    sequence_onehot = torch.zeros(4, self.max_size)
    sequence_onehot[:, :size] = one_hot_encode(sequence, "GAUC")

    mask = torch.zeros(self.max_size, dtype=torch.bool)
    mask[:size] = True
    bonds = torch.zeros(self.max_size, dtype=torch.long)
    causal_bonds = torch.zeros(self.max_size, dtype=torch.long)
    b, cb = block_bonds(structure)
    bonds[:size] = b
    causal_bonds[:size] = cb

    inputs = (
      sequence_onehot, bonds, causal_bonds, mask
    )
    outputs = (
      sequence_gt, mask
    )

    return inputs, outputs

class RNABlockRefoldRecomputeData(RNAGraphFoldRecomputeData):
  def __init__(self, min_size=10, max_size=200, pmut=0.1):
    self.min_size = min_size
    self.max_size = max_size
    self.pmut = pmut

  def __getitem__(self, index):
    size = random.randrange(self.min_size, self.max_size)
    seq = random.choices(["G", "A", "U", "C"], k=size)
    target_sequence = "".join(seq)
    step = random.randrange(1, 100)
    final_probability = 1 - (1 - self.pmut) ** (step - 1)
    step_sequence = "".join([
      random.choice("GAUC") if random.random() < final_probability else item
      for item in seq
    ])
    sequence = "".join([
      random.choice("GAUC") if random.random() < self.pmut else item
      for item in step_sequence
    ])
    target_structure, _ = fold(target_sequence)
    current_structure, _ = fold(sequence)

    sequence_gt = torch.zeros(self.max_size, dtype=torch.long)
    sequence_gt[:size] = one_hot_encode(step_sequence, "GAUC", numeric=True)
    sequence_onehot = torch.zeros(4, self.max_size)
    sequence_onehot[:, :size] = one_hot_encode(sequence, "GAUC")

    mask = torch.zeros(self.max_size, dtype=torch.bool)
    mask[:size] = True
    bonds = torch.zeros(self.max_size, dtype=torch.long)
    target = torch.zeros(self.max_size, dtype=torch.long)
    b, _ = block_bonds(current_structure)
    bonds[:size] = b
    t, _ = block_bonds(target_structure)
    target[:size] = t

    step = torch.tensor(step, dtype=torch.float)
    step = torch.tensor([
      (2 ** idx * step / 1000).sin()
      for idx in range(10)
    ], dtype=torch.float)

    inputs = (
      sequence_onehot, step, bonds, target, mask
    )
    outputs = (
      sequence_gt, mask
    )

    return inputs, outputs

class RNAGraphRefoldData(RefoldData):
  def __getitem__(self, index):
    source, source_structure, _, sequence, structure, _ = super().__getitem__(index)

    sequence_gt = one_hot_encode(sequence, "GAUC", numeric=True)
    source_onehot = one_hot_encode(source, "GAUC").permute(1, 0)
    sequence_onehot = one_hot_encode(sequence, "GAUC").permute(1, 0)
    source_bond, source_stem, source_unbound = tensor_structure(source_structure)
    bond, stem, unbound = tensor_structure(structure)
    rna = SubgraphStructure(torch.zeros(len(sequence), dtype=torch.long))

    inputs = (
      PackedTensor(source_onehot), rna,
      source_bond, bond,
      source_stem, stem,
      source_unbound, unbound
    )
    outputs = PackedTensor(sequence_gt, split=False)

    return inputs, outputs
