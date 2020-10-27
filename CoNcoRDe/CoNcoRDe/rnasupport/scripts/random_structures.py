from argparse import ArgumentParser
from rnasupport.samplers import fixed_length_structure

def parse_args():
  parser = ArgumentParser(
    description='Generate random RNA secondary structures of fixed length. ' \
                'Powered by the operad of RNA structures and its resource algebra.'
  )
  parser.add_argument('length', metavar='LENGTH', type=int,
                      help='Length of desired RNA structures.')
  parser.add_argument('--count', type=int, nargs='?', default=1,
                      help='Number of RNA structures to generate.')
  parser.add_argument('--min-E', type=int, nargs='?', default=1,
                      help='Minimum length of non-hairpin unpaired bases.')
  parser.add_argument('--max-E', type=int, nargs='?', default=10,
                      help='Maximum length of non-hairpin unpaired bases.')
  parser.add_argument('--min-H', type=int, nargs='?', default=3,
                      help='Minimum length of hairpin unpaired bases.')
  parser.add_argument('--max-H', type=int, nargs='?', default=10,
                      help='Maximum length of hairpin unpaired bases.')
  parser.add_argument('--min-S', type=int, nargs='?', default=2,
                      help='Minimum length of stems.')
  parser.add_argument('--max-S', type=int, nargs='?', default=10,
                      help='Maximum length of stems.')
  parser.add_argument('--min-Z', type=int, nargs='?', default=1,
                      help='Minimum length of continued stems.')
  parser.add_argument('--max-Z', type=int, nargs='?', default=10,
                      help='Maximum length of continued stems.')
  parser.add_argument('--min-M', type=int, nargs='?', default=0,
                      help='Minimum size of multiloops.')
  parser.add_argument('--max-M', type=int, nargs='?', default=5,
                      help='Maximum size of multiloops.')
  parser.add_argument('--min-B', type=int, nargs='?', default=1,
                      help='Minimum size of bulges.')
  parser.add_argument('--max-B', type=int, nargs='?', default=10,
                      help='Maximum size of bulges.')

  return parser.parse_args()


if __name__ == "__main__":
  opt = parse_args()
  cost = {
    "E": (opt.min_E, opt.max_E),
    "H": (opt.min_H, opt.max_H),
    "S": (opt.min_S, opt.max_S),
    "Z": (opt.min_Z, opt.max_Z),
    "M": (opt.min_M, opt.max_M),
    "B": (opt.min_B, opt.max_B),
  }

  for _ in range(opt.count):
    print(fixed_length_structure(opt.length, cost=cost))
