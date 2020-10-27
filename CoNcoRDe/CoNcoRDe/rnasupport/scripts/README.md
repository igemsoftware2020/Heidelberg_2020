# Scripts

Here, we provide a set of scripts for working with RNA secondary structures using methods from category theory.
To run these scripts, you will have to install `rnasupport` to your local python environment (`python setup.py develop`).

## RNA Secondary Structure Generation

We provide the script `random_structures.py` as a means to generate fixed-length RNA secondary structures adhering to user
constraints for minimum and maximum sizes of stems, loops, bulges and multiloops. After installing `rnasupport`, usage is simple:

```
usage: random_structures.py [-h] [--count [COUNT]] [--min-E [MIN_E]]
                            [--max-E [MAX_E]] [--min-H [MIN_H]]
                            [--max-H [MAX_H]] [--min-S [MIN_S]]
                            [--max-S [MAX_S]] [--min-Z [MIN_Z]]
                            [--max-Z [MAX_Z]] [--min-M [MIN_M]]
                            [--max-M [MAX_M]] [--min-B [MIN_B]]
                            [--max-B [MAX_B]]
                            LENGTH

Generate random RNA secondary structures of fixed length. Powered by the
operad of RNA structures and its resource algebra.

positional arguments:
  LENGTH           Length of desired RNA structures.

optional arguments:
  -h, --help       show this help message and exit
  --count [COUNT]  Number of RNA structures to generate.
  --min-E [MIN_E]  Minimum length of non-hairpin unpaired bases.
  --max-E [MAX_E]  Maximum length of non-hairpin unpaired bases.
  --min-H [MIN_H]  Minimum length of hairpin unpaired bases.
  --max-H [MAX_H]  Maximum length of hairpin unpaired bases.
  --min-S [MIN_S]  Minimum length of stems.
  --max-S [MAX_S]  Maximum length of stems.
  --min-Z [MIN_Z]  Minimum length of continued stems.
  --max-Z [MAX_Z]  Maximum length of continued stems.
  --min-M [MIN_M]  Minimum size of multiloops.
  --max-M [MAX_M]  Maximum size of multiloops.
  --min-B [MIN_B]  Minimum size of bulges.
  --max-B [MAX_B]  Maximum size of bulges.
```

For example, to generate 100 secondary structures using the standard size parameters, run the following:

```
python random_structures.py 100 --count 100
```

## Simple RNA Sequence Design

To design RNA sequences conditioned on a given RNA secondary structure and sequence constraints, we provide the general
`design_simple.py` script. This implements a simple genetic algorithm to search sequence space for RNA sequences conforming
both to sequence and structural constraings. Here too, usage is simple once `rnasupport` has been installed:

```
usage: design_simple.py [-h] [--size [SIZE]] [--fraction [FRACTION]]
                        [--random-fraction [RANDOM_FRACTION]]
                        [--max-steps [MAX_STEPS]] [--pmut [PMUT]]
                        SEQUENCE STRUCTURE

Search RNA sequences folding a given target structure.

positional arguments:
  SEQUENCE              RNA sequence constraints.
  STRUCTURE             Target RNA secondary structure.

optional arguments:
  -h, --help            show this help message and exit
  --size [SIZE]         size of the sequence pool.
  --fraction [FRACTION]
                        number of top sequences to propagate.
  --random-fraction [RANDOM_FRACTION]
                        number of sequences to randomly repopulate.
  --max-steps [MAX_STEPS]
                        maximum number of search steps.
  --pmut [PMUT]         probability of single-base mutation.
```

Just specify your target sequence and structure and let the algorithm run for some time. Should the algorithm fail to
find a satisfying solution, consider rerunning with larget `--size` and `--max-steps`.

### Flanking Hammer-Head Ribozyme (HHR) Design
As our project deals with the *in vivo* expression of functional RNA, a common design task features the design
of flanking Hammer-Head Ribozymes (HHR) for a given target RNA sequence. We have specifically provided a version
of our `design_simple.py` script specialized to this task. The script, `design_hhr.py` uses the same machinery
as `design_simple.py` to design flanking HHR for a sequence given just the 5' and 3' ends of that sequence.

```
usage: design_hhr.py [-h] [--size [SIZE]] [--fraction [FRACTION]]
                     [--random-fraction [RANDOM_FRACTION]]
                     [--max-steps [MAX_STEPS]] [--pmut [PMUT]]
                     5PRIME 3PRIME

Search flanking ribozymes for a target RNA sequence.

positional arguments:
  5PRIME                5' part of the target RNA (at least 6 bases).
  3PRIME                3' part of the target RNA (at least 4 bases).

optional arguments:
  -h, --help            show this help message and exit
  --size [SIZE]         size of the sequence pool.
  --fraction [FRACTION]
                        number of top sequences to propagate.
  --random-fraction [RANDOM_FRACTION]
                        number of sequences to randomly repopulate.
  --max-steps [MAX_STEPS]
                        maximum number of search steps.
  --pmut [PMUT]         probability of single-base mutation.
```

For example, you could use it like so:

```
python design_hhr.py GATTACA GUACAC
```

## CoNCoRDe

During our project, we have taken a view of RNA-structures grounded in category theory. As part of that approach,
we implemented our algorithm for Compositional Nested Conditional RNA Design (CoNCoRDe) as a more powerful alternative
to standard RNA sequence design algorithms. CoNCoRDe decomposes a target RNA structure into a partially ordered set
of substructures which it then attempts to recursively solve and re-compose into solutions for larger subproblems until
a solution to the target RNA structure is found. To this end, it stores substructure solutions in a library, which can
be reused accross RNA designs, resulting in a learning algorithm in the classical sense and accelerating the design
of families of target RNA structures.

To use CoNCoRDe from the command-line, we provided a script `concorde.py`, which is similar in usage to `design_simple.py`.

```
usage: concorde.py [-h] [--path [PATH]] [--size [SIZE]]
                   [--inner-size [INNER_SIZE]] [--fraction [FRACTION]]
                   [--random-fraction [RANDOM_FRACTION]]
                   [--max-steps [MAX_STEPS]] [--inner-steps [INNER_STEPS]]
                   [--regen [REGEN]] [--early]
                   STRUCTURE

Search RNA sequences folding a given target structure. Based on the operad of
RNA structures, decomposes the target structure into substructures which are
solved recursively.

positional arguments:
  STRUCTURE             Target RNA secondary structure.

optional arguments:
  -h, --help            show this help message and exit
  --path [PATH]         Path to save and load substructure library from.
  --size [SIZE]         size of the sequence pool.
  --inner-size [INNER_SIZE]
                        size of the sequence pool for substructure search.
  --fraction [FRACTION]
                        number of top sequences to propagate.
  --random-fraction [RANDOM_FRACTION]
                        number of sequences to randomly repopulate.
  --max-steps [MAX_STEPS]
                        maximum number of search steps.
  --inner-steps [INNER_STEPS]
                        maximum number of search steps for substructure
                        search.
  --regen [REGEN]       probability of diversifying the substructure library
                        at each step.
  --early               Finish sequence search early, if matching sequence is
                        found?
```

You can provide a custom path to save your CoNCoRDe substructure library for subsequent runs and limit computational
resources in the inner loop of substructure search.

To design, for example a difficult but small tetraloop, you could run the following command:

```
python rnasupport/scripts/concorde.py "((....))((....))((....))((....))" --max-steps 1000
```

If the first attempt at design fails, the substructure library makes restarting a run trivially easy.
Just run the same command again, and CoNCoRDe will load all previously solved substructures first and continue sequence search.

## Design RNA Triple-Helix

During our project, we have taken a view of RNA-structures grounded in category theory. As part of that approach,
we implemented our algorithm for Compositional Nested Conditional RNA Design (CoNCoRDe) as a more powerful alternative
to standard RNA sequence design algorithms. CoNCoRDe decomposes a target RNA structure into a partially ordered set
of substructures which it then attempts to recursively solve and re-compose into solutions for larger subproblems until
a solution to the target RNA structure is found. To this end, it stores substructure solutions in a library, which can
be reused accross RNA designs, resulting in a learning algorithm in the classical sense and accelerating the design
of families of target RNA structures.

To use our Design tool for RNAs in a Triple-Helix from the command-line, we provided a script `design_rna_triple_helix.py`:

```
usage: design_rna_triple_helix.py [-h] [--path [PATH]] [--size [SIZE]]
                   [--dna_sequence [DNA_SEQUENCE]]

Translate DNA in Triple-Helix RNA based on Table 1 in Article: Kunkler, Charlotte N and Hulewicz, 
Jacob P and Hickman, Sarah C and Wang, Matthew C and McCown, Phillip J and Brown, Jessica A (2019) 
{Stability of an RNA•DNA–DNA triple helix depends on base triplet composition and length of the RNA 
third strand}.Nucleic Acids Research 47, 7213-7222.')

optional arguments:
  -h, --help            show this help message and exit

```

To generate an RNA to use for RNA Triple-Helix design, you could run for example:

```
python rnasupport/scripts/design_rna_triple_helix.py --dna_sequence "AGGT"
```