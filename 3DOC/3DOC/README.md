# 3DOC

## Overview

We are introducing a pipeline to create protein database format (PDB) files out of concatenated protein sequences, which can be interconnected via amino acid linkers. The pipeline enables creating fusion protein-binding domains and fusion proteins. For instance, it can be used to create Pumby and PPR protein sequences and PDBs, which can be linked modularly to bind specific RNAs. ((Adamala AP *et al.*, 2016) (Coquille S *et al.*, 2014)). The PDB file creation pipeline uses BLASTp (https://blast.ncbi.nlm.nih.gov/Blast.cgi) together with PyRosetta respectively trRosetta to generate a PDB files for protein sequences ((Chaudhury S *et al.*, 2010); (Yang J *et al.*, 2020)).  If several sequences are inputted, the PDBs are fused or interconnected automatically with the inputted Amino Acid Linker via PyRosetta. We are also offering two scripts "rna_denovo_preparation.py" respectively "rnp_structure_prediction_preparation.py", which prepare the output of 3DOC for RNA denovo (FARFAR) and RNP Structure prediction protocols in Rosetta ((Das R *et al.*, 2010), (Cheng C *et al.*, 2015), (Kappel K and Das R, 2019)).

The full documentation can be seen accessed on the [iGEM Heidelberg 2020 Wiki](https://2020.igem.org/Team:Heidelberg/Software/3DOC).

## Installation:

Run setup.py in the root directory like this:

```
python setup.py devel
```

## System Requirements:
* Python 3.7
*	Bash Version 4+
*	PyRosetta-3.7.Release (Tested on Scientific/Red Hat Linux Version.)

## Packages:
*	tensorflow 1.13+
*	biopython 1.78+
* ViennaRNA

## Bibliography:

[1] Adamala, K.P., Martin-Alarcon, D., & Boyden, E. (2016). Programmable RNA-binding protein composed of repeats of a single modular unit. Proceedings of the National Academy of Sciences, 113, E2579 - E2588.

[2] Chaudhury, S., Lyskov, S., & Gray, J.J. (2010). PyRosetta: a script-based interface for implementing molecular modeling algorithms using Rosetta. Bioinformatics, 26 5, 689-91.

[3] Cheng, C., Chou, F., & Das, R. (2015). Modeling complex RNA tertiary folds with Rosetta. Methods in enzymology, 553, 35-64.

[4] Coquille, S.C., Filipovska, A., Chia, T., Rajappa, L., Lingford, J.P., Razif, M., Thore, S., & Rackham, O. (2014). An artificial PPR scaffold for programmable RNA recognition. Nature communications, 5, 5729.

[5] Das, R., Karanicolas, J., & Baker, D. (2010). Atomic accuracy in predicting and designing non-canonical RNA structure. Nature methods, 7, 291 - 294.

[6] Kappel, K., & Das, R. (2019). Sampling Native-like Structures of RNA-Protein Complexes through Rosetta Folding and Docking. Structure, 27 1, 140-151.e5.
