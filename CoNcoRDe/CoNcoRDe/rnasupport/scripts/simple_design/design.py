import random
import RNA

CONSTRAINTS = dict(
  N="GAUC",
  Y="UC",
  R="GA",
  H="AUC",
  B="GUC",
  V="GAC",
  D="GAU",
  G="G",
  A="A",
  U="U",
  C="C"
)

COMPLEMENT = dict(
  U="GA",
  A="U",
  G="UC",
  C="G"
)

def parse_structure(structure):
  """Parses an RNA secondary structure into a
  dictionary format.

  Args:
      structure (str): dot-bracket structure string.

  Returns:
      Dictionary containing pairs of paired base indices.
  """
  stack = []
  result = {}
  for idx, val in enumerate(structure):
    if val == "(":
      stack.append(idx)
    if val == ")":
      idy = stack.pop()
      result[idx] = idy
      result[idy] = idx
  return result

def _restrict_complement(base, constraint):
  return [
    comp
    for comp in COMPLEMENT[base]
    if comp in constraint
  ]

def _pair_constraints(constraint, comp_constraint):
  result = []
  for cc in constraint:
    for cc_p in COMPLEMENT[cc]:
      if cc_p in comp_constraint:
        result.append(cc)
  return result
              
def random_sequence(constraint, structure):
  """Generates a random RNA sequence roughly corresponding
  to an extended base alphabet constraint string coupled with
  a target secondary structure.

  Args:
    constraint (str): extended base alphabet RNA sequence.
      Constrains the generated sequence to adhere to the
      sequence specification given by this string of extended
      bases.
    structure (str): target dot-bracket structure.
  """
  pairs = parse_structure(structure)
  result = ""
  for idx, cc in enumerate(constraint):
    if idx in pairs and pairs[idx] < idx:
      comp = _restrict_complement(
        result[pairs[idx]],
        CONSTRAINTS[constraint[idx]]
      )
      result += random.choice(comp)
    elif idx in pairs:
      pair_cstr = CONSTRAINTS[constraint[pairs[idx]]]
      cstr = CONSTRAINTS[constraint[idx]]
      comp = _pair_constraints(
        cstr,
        pair_cstr
      )
      result += random.choice(comp)
    else:
      result += random.choice(CONSTRAINTS[cc])
  return result

def mutate_sequence(sequence, constraint, structure, pmut=0.1):
  """Mutates a given RNA sequence roughly adhering to
  an extended base alphabet constraint string coupled with
  a target secondary structure.

  Args:
    constraint (str): extended base alphabet RNA sequence.
      Constrains the generated sequence to adhere to the
      sequence specification given by this string of extended
      bases.
    structure (str): target dot-bracket structure.
    pmut (float): probability of single-base mutation.
  """
  pairs = parse_structure(structure)
  rpos = random.choice(range(len(sequence)))
  result = ""
  for idx, cc in enumerate(constraint):
    if pmut < random.random():
      result += sequence[idx]
    elif idx in pairs and pairs[idx] < idx:
      comp = _restrict_complement(
        result[pairs[idx]],
        CONSTRAINTS[constraint[idx]]
      )
      result += random.choice(comp)
    elif idx in pairs:
      pair_cstr = CONSTRAINTS[constraint[pairs[idx]]]
      cstr = CONSTRAINTS[constraint[idx]]
      comp = _pair_constraints(
        cstr,
        pair_cstr
      )
      result += random.choice(comp)
    else:
      result += random.choice(CONSTRAINTS[cc])
  return result

def loop(length):
  """Constructs the dot-bracket structure of a loop.

  Args:
    length (int): length of the loop.
  """
  return "." * length

def stem(length, structure):
  """Constructs the dot-bracket structure of a stem.

  Args:
    length (int): length of the stem.
    structure (str): dot-bracket structure flanked by the stem.
  """
  return "(" * length + structure + ")" * length

def _init_pool(sequence, structure, size=20):
  return [
    random_sequence(sequence, structure)
    for _ in range(size)
  ]

def _update_pool(pool, fitness, sequence, structure,
                 fraction=4, random_fraction=5, pmut=0.1):
  fit = list(map(fitness, pool))
  keys = sorted(list(range(len(fit))), key=lambda x: fit[x])
  seed = [
    pool[keys[idx]]
    for idx in range(fraction)
  ]
  new_pool = [
    mutate_sequence(random.choice(seed), sequence, structure, pmut=pmut)
    for _ in range(len(pool) - random_fraction)
  ] + _init_pool(sequence, structure, size=random_fraction)
  return new_pool

def _structure_fitness(structure):
  def _target(sequence):
    predicted, energy = RNA.fold(sequence)
    mismatches = sum([
      int(x != y)
      for x, y in zip(predicted, structure)
    ])
    return 100 * mismatches + energy
  return _target

def design(sequence, structure, size=20, fraction=4,
           random_fraction=5, max_steps=100, pmut=0.1):
  """Searches RNAs folding a given structure and adhering to sequence constraints
  using a simple genetic algorithm.

  Args:
    sequence (str): expanded alphabet constraint sequence.
    structure (str): dot-bracket target secondary structure.
    size (int): size of the sequence pool.
    fraction (int): number of top sequences to use for pool
      repopulation.
    random_fraction (int): number of sequences to repopulate
      uniformly at random.
    max_steps (int): maximum number of search steps.
    pmut (float): single-base mutation probability.
  """
  fitness = _structure_fitness(structure)
  pool = _init_pool(sequence, structure, size=size)
  best = None
  best_value = None
  for idx in range(max_steps):
    values = list(map(fitness, pool))
    local_best, local_value = sorted(list(zip(pool, values)), key=lambda x: x[1])[0]
    if not best or local_value < best_value:
      best = local_best
      best_value = local_value
      print(best, best_value)
    pool = _update_pool(
      pool, fitness, sequence,
      structure, fraction=fraction,
      random_fraction=random_fraction
    )
  return best, best_value
