# Preprocessing Pipeline for Geneious Prime/ Manual for 3DOC

We are presenting the Geneious Prime workflow “3DOC” as well as instructions for a manual pipeline for preparation of DNA sequence files for the PDB file creation pipeline. Also, we are offering a library of Pumby and PPR modules as well as common RNA linkers for optional usage in the Geneious Prime workflow respectively the manual pipeline. The files can be imported into Geneious Prime via “File (In the menu.) -> Import -> From File”.

## Preparation for PDB file creation pipeline

The Geneious Prime workflow is offered in the file "workflow_3DOC_preprocessing.geneiousWorkflow" and needs to be imported into Geneious Prime via “Tools (In the menu.) -> Workflows -> Manage Workflows -> Import”. In Geneious Prime the DNA-sequence files need to be prepared for the workflow. The preparation steps as well as some further steps apply for the manual preparation of sequence files as well:

1. Please use the folder “3DOC” for Geneious Prime or as ZIP file for manual preparation offered in the Github. Copy all DNA-sequence files, which you want to work with in this pipeline, into the folder "3DOC". YOu may also copy pumby (pumilio) and ppr protein modules from their corresponding folders.

1. Change the names of DNA-sequence documents, so they are all written in the format "X sequence_name". The words are interconnected with an underscore.

1. The sequences are numbered with ascending numbers from the beginning to the last sequence in the RBP-binding fusion domains respectively fusion protein. The first protein is marked by adding "- begin" at the end of the dna-sequence document name. See the example:

		"1 protein_2 - begin"
		"2 protein_1"
		"3 protein_3"

1. Protein linker can be copied from the folder "Linker" or inserted into Geneious on your own. They are numbered according to the convention below, to demonstrate, which sequences they interconnect. The word "linker" in the file is necessary. The smaller protein sequence number is written at the front with a dash connecting the smaller protein sequence number. Only protein linkers between proteins, which are direct neighbours, are allowed. See the example:

		"1-2 linker_1"

1. Removal of Methionine at the begin of sequences.
markdown how to

  	1. Geneious Prime: The folder "Trimming" is necessary to perform the trimming of the DNA-sequences for Methione at the beginning of the sequences except for the first sequence. We cannot guarantee for correct functioning if the folder is modified.
  	1. Manual preparation: Please delete the bases “ATG” (Are translated to Methionine.) for all sequences except for the beginning sequence, if it can be found at 5’ end.

1. Please delete all STOP codons at the end of all sequences except for the last sequence. They may occur for Sequences from iGEM Registry and SynBioHub.

1. Output
	1. Geneious Prime: The workflow outputs the translation of the sequences directly to a folder, which can be chosen by the user. It also outputs the concatenated DNA sequence.
	1. Manual preparation: You may translate the sequences separately with another tool (e.g. EMBOSS by EMBL [URL: https://www.ebi.ac.uk/Tools/st/emboss_transeq/])  and concatenate the DNA-sequences manually.

## Library of Pumby (Pumilio) and PPR Modules

Both pumby (pumilio) and ppr modules in the folders "Pumby (Pumilio) Modules" as well as "PPR Modules" can be used for the Geneious workflow as well as the manual pipeline. Both types of protein modules can be combined in a modular way and bind specific to a certain nucleotide. 

The pumby (Pumilio) modules are based on a consensus sequence of Pumilio family. When using the pumby modules please refer to the article “Programmable RNA-binding protein composed of repeats of a single modular unit”, published in 2016 in the journal PNAS (Adamala et al., 2016). Please note, that in accordance to the authors a sequence of pumby modules is supposed to start with the unit “X pumby_module_start_NUCLEOTIDE” and to end with the unit “X pumby_module_end_NUCLEOTIDE”. In accordance to experiments presented in the paper the pipeline will inform the user, if the “start” and “end” module are positioned in the wrong position and if the amount of modules is wrong (Sequence Length: 6/ 10/ 12/ 18 or 24)

The PPR Modules are split into two groups:

1. “X cPPR”-polyNUCLEOTIDE” group, which is based on a consensus PPR protein family motif. When using this group of PPR modules please refer to the article “An artificial PPR scaffold for programmable RNA recognition”, published in 2014 in the journal Nature Communications (Coquille et al., 2014). 
1. “X ppr10_repeat6/7_NUCLEOTIDE” group, which is based on single amino acid substitutions of the repeat 6 respectively 7 of the PPR10 protein. When using this group of PPR modules please refer to the article “A combinatorial amino acid code for RNA recognition by pentatricopeptide repeat proteins.”, published in 2012 in the journal PLoS Genetics (Barkan et al., 2012). 

We also provide Amino Acid Linkers described in literature in the folder "Linker" to interconnect protein domains or fusion protein parts with each other (Chen X et al., 2013).
