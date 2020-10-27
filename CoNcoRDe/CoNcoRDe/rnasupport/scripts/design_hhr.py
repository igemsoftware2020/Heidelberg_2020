#!/usr/bin/env python

from argparse import ArgumentParser

from rnasupport.scripts.simple_design.design import stem, loop, design

def _hammerhead_template_5_structure(stem_1=6, stem_2=3, stem_3=4):
  stem_2 = stem(stem_2 + 1, loop(4))
  stem_3 = stem(stem_3 + 1, loop(4))
  structure = loop(7) + stem_2 + loop(3) + stem_3 + loop(1)
  structure = stem(stem_1, structure)
  return structure

def _hammerhead_template_5_sequence(target, stem_2=3, stem_3=4):
  seq = "N" * len(target) + "CUGANGAR"
  seq += "N" * (stem_2 * 2 + 4) + "YGAAA"
  seq += "N" * (stem_3 * 2 + 4) + "UH"
  return seq + target

def hh_5_template(target, stem_2=3, stem_3=4):
  structure = _hammerhead_template_5_structure(
    stem_1=len(target), stem_2=stem_2, stem_3=stem_3
  )
  sequence = _hammerhead_template_5_sequence(
    target, stem_2, stem_3
  )
  return sequence, structure

def _hammerhead_template_3_structure(stem_1=6, stem_2=3, stem_3=4):
  stem_2 = stem(stem_2, loop(4))
  stem_3 = stem(stem_3 + 1, loop(4))
  structure = loop(1) + stem_2 + loop(7) + stem_3 + loop(3)
  structure = stem(stem_1 + 1, structure)
  return structure

def _hammerhead_template_3_sequence(target, stem_2=3, stem_3=4):
  seq = target + "UH"
  seq += "N" * (stem_2 * 2 + 4) + "CUGANGAR"
  seq += "N" * (stem_3 * 2 + 4) + "YGAAA" + "N" * len(target)
  return seq

def hh_3_template(target, stem_2=4, stem_3=3):
  structure = _hammerhead_template_3_structure(
    stem_1=len(target), stem_2=stem_2, stem_3=stem_3
  )
  sequence = _hammerhead_template_3_sequence(
    target, stem_2, stem_3
  )
  return sequence, structure

def parse_args():
  parser = ArgumentParser(
    description='Search flanking ribozymes for a target RNA sequence.'
  )
  parser.add_argument('p5', metavar='5PRIME', type=str,
                      help='5\' part of the target RNA (at least 6 bases).')
  parser.add_argument('p3', metavar='3PRIME', type=str,
                      help='3\' part of the target RNA (at least 4 bases).')
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

  print("Searching 5' ribozymes...")
  result_5p = design(
    *hh_5_template(opt.p5.upper().replace("T", "U")),
    size=opt.size, fraction=opt.fraction,
    random_fraction=opt.random_fraction,
    max_steps=opt.max_steps, pmut=opt.pmut
  )

  print("Searching 3' ribozymes...")
  result_3p = design(
    *hh_3_template(opt.p3.upper().replace("T", "U")),
    size=opt.size, fraction=opt.fraction,
    random_fraction=opt.random_fraction,
    max_steps=opt.max_steps, pmut=opt.pmut
  )

  print("Results:")
  print("5' HHR:", result_5p[0], "Energy:", result_5p[1])
  print("3' HHR:", result_3p[0], "Energy:", result_3p[1])
