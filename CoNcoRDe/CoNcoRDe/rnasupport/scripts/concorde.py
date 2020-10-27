#!/usr/bin/env python

import os
from argparse import ArgumentParser

from rnasupport.hierarchy import HierarchyDesign, parse_hierarchy

def parse_args():
  parser = ArgumentParser(
    description='Search RNA sequences folding a given target structure. ' \
                'Based on the operad of RNA structures, decomposes the target ' \
                'structure into substructures which are solved recursively.'
  )
  parser.add_argument('structure', metavar='STRUCTURE', type=str,
                      help='Target RNA secondary structure.')
  parser.add_argument('--path', type=str, nargs='?', default=".concorde-library.dat",
                      help='Path to save and load substructure library from.')
  parser.add_argument('--size', type=int, nargs='?', default=100,
                      help='size of the sequence pool.')
  parser.add_argument('--inner-size', type=int, nargs='?', default=10,
                      help='size of the sequence pool for substructure search.')
  parser.add_argument('--fraction', type=int, nargs='?', default=5,
                      help='number of top sequences to propagate.')
  parser.add_argument('--random-fraction', type=int, nargs='?', default=5,
                      help='number of sequences to randomly repopulate.')
  parser.add_argument('--max-steps', type=int, nargs='?', default=500,
                      help='maximum number of search steps.')
  parser.add_argument('--inner-steps', type=int, nargs='?', default=5,
                      help='maximum number of search steps for substructure search.')
  parser.add_argument('--regen', type=float, nargs='?', default=0.01,
                      help='probability of diversifying the substructure library at each step.')
  parser.add_argument('--early', dest='early', action='store_true',
                      help='Finish sequence search early, if matching sequence is found?')

  return parser.parse_args()

if __name__ == "__main__":
  opt = parse_args()

  design = HierarchyDesign(
    population_size=opt.size,
    library_regen=opt.regen,
    steps=opt.max_steps,
    inner_steps=opt.inner_steps,
    inner_population_size=opt.inner_size
  )

  if os.path.isfile(opt.path):
    design.load_library(opt.path)

  try:
    print("Searching RNAs...")
    result = design.design(
      parse_hierarchy(opt.structure),
      steps=opt.max_steps,
      population_size=opt.size,
      early=opt.early
    )

    print("RNA:", result[0], "Structure:", result[1])
  except Exception as e:
    raise e
  finally:
    design.dump_library(opt.path)
