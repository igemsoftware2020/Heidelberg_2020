import torch

from torchsupport.structured import ScatterStructure

def _chase_cycle(idx, bond_dict, size):
  current_list = []
  start = idx
  current = idx + 1
  jumped = False
  if idx in bond_dict:
    current = bond_dict[idx]
    jumped = True
  current_list += [start]
  while current != start and current < size:
    current_list.append(current)
    if current in bond_dict and not jumped:
      current = bond_dict[current]
      jumped = True
    else:
      current = current + 1
      jumped = False
  return tuple(sorted(list(set(current_list))))

def parse_secondary_structure(structure):
  size = len(structure)
  bond_dict = {}
  unbound_contiguous_list = []
  stem_contiguous_list = []
  stack = []
  for idx, char in enumerate(structure):
    if char == "(":
      stack.append(idx)
    if char == ")":
      partner = stack.pop()
      bond_dict[idx] = partner
      bond_dict[partner] = idx

  unbound = False
  current_stem = []
  touched = []
  for idx in range(size):
    if unbound and idx in bond_dict:
      unbound = False
      unbound_contiguous_list.append(_chase_cycle(idx, bond_dict, size))
    elif not unbound and idx not in bond_dict:
      unbound = True
      if current_stem:
        stem_contiguous_list.append(tuple(list(set(current_stem))))
      current_stem = []
      unbound_contiguous_list.append(_chase_cycle(idx, bond_dict, size))
    if not unbound and idx in bond_dict and idx not in touched:
      current_stem.append(idx)
      current_stem.append(bond_dict[idx])
      touched.append(idx)
      touched.append(bond_dict[idx])
  return bond_dict, list(set(stem_contiguous_list)), list(set(unbound_contiguous_list))

def _bond_tensors(bond_dict, size):
  bond_neighbours = []
  bond_indices = []

  for idx in range(size):
    if idx in bond_dict:
      bond_indices += [idx, idx]
      bond_neighbours += [idx, bond_dict[idx]]

  bond_neighbours = torch.tensor(bond_neighbours, dtype=torch.long)
  bond_indices = torch.tensor(bond_indices, dtype=torch.long)

  return bond_neighbours, bond_indices

def _stem_tensors(stems, size):
  stem_neighbours = []
  stem_indices = []
  for stem in stems:
    for index in stem:
      for other_index in stem:
        stem_indices.append(index)
        stem_neighbours.append(other_index)

  indices = list(range(len(stem_indices)))
  indices = sorted(indices, key=lambda x: stem_indices[x])
  stem_neighbours = [stem_neighbours[index] for index in indices]
  stem_indices = [stem_indices[index] for index in indices]

  stem_neighbours = torch.tensor(stem_neighbours, dtype=torch.long)
  stem_indices = torch.tensor(stem_indices, dtype=torch.long)

  return stem_neighbours, stem_indices

def tensor_structure(structure):
  size = len(structure)
  bond_dict, stems, unbounds = parse_secondary_structure(structure)

  bond_tensors = _bond_tensors(bond_dict, size)
  stem_tensors = _stem_tensors(stems, size)
  unbound_tensors = _stem_tensors(unbounds, size)

  bond_structure = ScatterStructure(0, 0, *bond_tensors, node_count=size)
  stem_structure = ScatterStructure(0, 0, *stem_tensors, node_count=size)
  unbound_structure = ScatterStructure(0, 0, *unbound_tensors, node_count=size)

  return bond_structure, stem_structure, unbound_structure

def block_bonds(structure):
  size = len(structure)
  bond_dict, stems, unbounds = parse_secondary_structure(structure)
  result = torch.zeros(size, dtype=torch.long)
  causal_result = torch.zeros(size, dtype=torch.long)
  for key in range(size):
    if key in bond_dict:
      value = bond_dict[key]
      result[key] = value
      if value < key:
        causal_result[key] = value
      else:
        causal_result[key] = key
    else:
      result[key] = key
      causal_result[key] = key
  return result, causal_result
