#!/usr/bin/env python

from argparse import ArgumentParser

from rnasupport.scripts.simple_design.design import design

def parse_args():
  parser = ArgumentParser(
    description='Search RNA sequences folding a given target structure.'
  )
  parser.add_argument('sequence', metavar='SEQUENCE', type=str,
                      help='RNA sequence constraints.')
  parser.add_argument('structure', metavar='STRUCTURE', type=str,
                      help='Target RNA secondary structure.')
  parser.add_argument('--size', type=int, nargs='?', default=100,
                      help='size of the sequence pool.')
  parser.add_argument('--fraction', type=int, nargs='?', default=5,
                      help='number of top sequences to propagate.')
  parser.add_argument('--random-fraction', type=int, nargs='?', default=5,
                      help='number of sequences to randomly repopulate.')
  parser.add_argument('--max-steps', type=int, nargs='?', default=20,
                      help='maximum number of search steps.')
  parser.add_argument('--pmut', type=float, nargs='?', default=0.1,
                      help='probability of single-base mutation.')

  return parser.parse_args()

if __name__ == "__main__":
  opt = parse_args()

  print("Searching RNAs...")
  result = design(
    opt.sequence.upper().replace("T", "U"), opt.structure,
    size=opt.size, fraction=opt.fraction,
    random_fraction=opt.random_fraction,
    max_steps=opt.max_steps, pmut=opt.pmut
  )

  print("Results:")
  print("RNA:", result[0], "Energy:", result[1])
