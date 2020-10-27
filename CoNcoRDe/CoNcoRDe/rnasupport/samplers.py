r"""This module defines a uniform random sampler of plausible RNA secondary structures
  based on the free operad of RNA secondary structures on the following set of generators:
    - Stems: $\{1, 2, 3, ...\} \to S: 1$
    - Loops: $\{3, 4, ...\} \to H: 0$
    - Multiloops: $\{2, 3, ...\} \to M: n$
    - Bulges: $\{0, 1, ...\} \times \{0, 1, ...\} \to B: 1$
"""

import random

class RNANode:
  r"""Operation in an operad of RNA structures. We chose the name 'Node' to allude to
    the tree-like structure of RNA secondary structure without pseudoknots.

  Args:
    kind ("E", "H", "S", "Z", "B", "M"): generator kind of the most recently composed operation.
    shape (int or tuple): shape of the most recently composed generator (e.g. stem length,
      loop size, bulge sizes, multiloop size).
    children (List[RNANode]): composed operations, i.e. children in the composition tree.
  """
  def __init__(self, kind, shape=None, children=None):
    self.kind = kind
    self.shape = shape
    self.children = children or []

  def stringify(self, depth=0):
    r"""Returns a string representation of a given RNANode."""
    indent = "  " * depth
    result = f"{indent}{self.kind}: {self.shape}\n"
    for child in self.children:
      result += child.stringify(depth=depth + 1)
    return result

  def as_tuple(self):
    r"""Returns a tuple representation of a given RNANode."""
    return (
      self.kind, self.shape, tuple(self.children)
    )

  def __hash__(self):
    return hash(self.as_tuple())

  def __repr__(self):
    return self.stringify()

def cost_bounds(kind, resource, cost=None):
  r"""Implementation of the resource coalgebra on the generators.
  Args:
    kind: generator kind.
    resource: available resource (i.e. sequence length).

  Returns:
    Minimum and maximum sequence length consumption for a given
    generator kind.
  """
  if kind == "E":
    return cost["E"][0], min(resource, cost["E"][1])
  if kind == "H":
    return resource, resource
  if kind in ["S", "Z"]:
    minimum = 2 * cost[kind][0] + cost["H"][0]
    maximum = min(cost[kind][1], (resource - cost["H"][0]) // 2)
    return minimum, maximum
  if kind == "B":
    stem_minimum = cost_bounds("Z", 0, cost=cost)[0]
    minimum = cost["B"][0] + stem_minimum
    maximum = min(cost["B"][1], (resource - stem_minimum))
    return minimum, maximum
  if kind == "M":
    stem_minimum = cost_bounds("S", 0, cost=cost)[0]
    minimum = 3 * cost["M"][0] + 2 * stem_minimum
    max_count = (resource - cost["M"][0]) // (stem_minimum + cost["M"][0])
    return minimum, max_count

def fixed_length_bulge(resource, cost=None):
  r"""Uses the resource coalgebra to calculate cost bounds for
  a random bulge structure and samples a plausible structure
  containing a bulge uniformly at random.

  Args:
    resource (int): maximum sequence length to invest.
    cost (Optional[dict]): dictionary of minimum and maximum
      resource consumption for all generators.

  Returns:
    Random RNANode representing an RNA secondary structure
    starting with a bulge.
  """
  min_bound = cost["B"][0]
  _, max_bound = cost_bounds("B", resource, cost=cost)

  left = random.randint(min_bound, max_bound)
  right = random.randint(0, max_bound - left)
  resource -= left + right
  shape = (left, right)

  result = RNANode("B", shape=shape)

  result.children.append(
    fixed_length_stem(resource, stem="Z", cost=cost)
  )

  return result

def fixed_length_junction(resource, cost=None):
  r"""Uses the resource coalgebra to calculate cost bounds for
  a random multiloop structure and samples a plausible structure
  containing a multiloop uniformly at random.

  Args:
    resource (int): maximum sequence length to invest.
    cost (Optional[dict]): dictionary of minimum and maximum
      resource consumption for all generators.

  Returns:
    Random RNANode representing an RNA secondary structure
    starting with a multiloop.
  """
  min_bound = cost["M"][0]
  _, max_count = cost_bounds("M", resource, cost=cost)
  count = min(5, random.randint(2, max_count))

  stem_minimum = cost_bounds("S", 0, cost=cost)[0]
  stem_resource_max = resource - (count + 1) * min_bound
  stem_resource_min = max(count * stem_minimum, resource - (count + 1) * cost["M"][1])

  stem_resource = random.randint(stem_resource_min, stem_resource_max)

  stems = []
  stem_consumption = []
  for idx in range(count):
    _, consumed = sample_allowed(
      stem_resource - (count - idx - 1) * stem_minimum,
      symbols=["S"],
      cost=cost
    )
    if idx == count - 1:
      consumed = stem_resource
    stem_consumption.append(consumed)
    resource -= consumed
    stem_resource -= consumed
    stems.append(
      fixed_length_stem(consumed, cost=cost)
    )
  shape = []
  assert stem_resource == 0
  for idx in range(count):
    consumed = random.randint(
      max(min_bound, resource - (count - idx) * cost["M"][1]),
      min(cost["M"][1], resource - (count - idx) * min_bound)
    )
    resource -= consumed
    shape.append(consumed)
  shape.append(resource)
  assert all(map(lambda x: x <= cost["M"][1], shape))
  resource = 0
  result = RNANode("M", shape=shape, children=stems)
  return result

def fixed_length_stem(resource, stem="S", cost=None):
  r"""Uses the resource coalgebra to calculate cost bounds for
  a random stem structure and samples a plausible structure
  containing a stem uniformly at random.

  Args:
    resource (int): maximum sequence length to invest.
    cost (Optional[dict]): dictionary of minimum and maximum
      resource consumption for all generators.

  Returns:
    Random RNANode representing an RNA secondary structure
    starting with a stem.
  """
  min_bound = cost[stem][0]
  _, max_bound = cost_bounds(stem, resource, cost=cost)

  consumed = random.randint(min_bound, max_bound)
  resource -= 2 * consumed

  result = RNANode(stem, shape=consumed)
  kind, consumed = sample_allowed(resource, symbols=["B", "M", "H"], cost=cost)

  if kind == "H":
    result.children.append(
      RNANode("H", shape=consumed)
    )
  elif kind == "M":
    result.children.append(
      fixed_length_junction(resource, cost=cost)
    )
  elif kind == "B":
    result.children.append(
      fixed_length_bulge(resource, cost=cost)
    )
  return result

def allowed_symbols(resource, symbols, cost):
  r"""Uses the resource coalgebra to find generators with satisfiable
  resource bounds for random generator selection.

  Args:
    resource (int): maximum sequence length to invest.
    symbols (list): list of possible generators to choose from.
    cost (Optional[dict]): dictionary of minimum and maximum
      resource consumption for all generators.

  Returns:
    List of generators with resource bounds less than or equal to the
    available resources.
  """
  allowed = []
  for symbol in symbols:
    min_bound, _ = cost_bounds(symbol, resource, cost=cost)
    if resource >= min_bound:
      bulge_minimum, _ = cost_bounds("B", resource, cost=cost)
      if symbol == "H" and resource >= bulge_minimum and resource >= cost["H"][1]:
        continue
      allowed.append(symbol)
  return allowed

def sample_allowed(resource, symbols=None, cost=None):
  r"""Uses the resource coalgebra to sample and compose a random
  generator with satisfiable resource bounds.

  Args:
    resource (int): maximum sequence length to invest.
    symbols (list): list of possible generators to choose from.
    cost (Optional[dict]): dictionary of minimum and maximum
      resource consumption for all generators.

  Returns:
    - pick: randomly chose admissible generator.
    - consumed: amount of consumed resources.
  """

  allowed = allowed_symbols(resource, symbols, cost)
  pick = random.choice(allowed)

  min_bound, _ = cost_bounds(pick, resource, cost=cost)
  max_bound = min(cost["E"][1], resource) if pick == "E" else resource
  consumed = random.randint(min_bound, max_bound)

  return pick, consumed

def fixed_length_map(resource, cost=None):
  r"""Uses the resource coalgebra on the operad of RNA secondary structures
  to sample an RNA secondary structure uniformly at random given a set of
  length-constraints on the generators.

  Args:
    resource (int): total desired sequence length.
    cost (Optional[dict]): dictionary of length constraints on the generators.
  """
  if cost is None:
    cost = {
      "E": (1, 10),
      "H": (3, 10),
      "S": (2, 10),
      "Z": (1, 10),
      "M": (0, 5),
      "B": (1, 10),
    }

  result = RNANode("R")
  while resource > 0:
    kind, consumed = sample_allowed(resource, symbols=["S"] * 4 + ["E"], cost=cost)
    resource -= consumed
    if kind == "E":
      result.children.append(
        RNANode("E", shape=consumed)
      )
    elif kind == "S":
      result.children.append(
        fixed_length_stem(consumed, cost=cost)
      )
  return result

def to_dot_bracket(tree):
  r"""Pretty-prints an RNANode as a dot-bracket string."""
  if tree.kind == 'R':
    result = "".join([
      to_dot_bracket(sub)
      for sub in tree.children
    ])
  elif tree.kind == 'E':
    result = "." * tree.shape
  elif tree.kind == 'H':
    result = "." * tree.shape
  elif tree.kind in ['S', 'Z']:
    inner = to_dot_bracket(
      tree.children[0]
    )
    count = tree.shape
    result = "(" * count + inner + ")" * count
  elif tree.kind == 'B':
    inner = to_dot_bracket(
      tree.children[0]
    )
    result = "." * tree.shape[0] + inner + "." * tree.shape[1]
  elif tree.kind == "M":
    subs = [
      to_dot_bracket(sub)
      for sub in tree.children
    ]
    spacers = [
      "." * size
      for size in tree.shape
    ]
    combined = [spacers[0]]
    for sub, spacer in zip(subs, spacers[1:]):
      combined.append(sub)
      combined.append(spacer)
    result = "".join(combined)
  return result

def fixed_length_structure(resource, cost=None):
  return to_dot_bracket(fixed_length_map(resource, cost=cost))
