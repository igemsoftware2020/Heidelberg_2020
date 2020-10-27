import os
import sys
sys.path.append(os.getcwd())
sys.path.append("./data")
sys.path.append("./tools")
sys.path.append("./Databases")
sys.path.append("./Databases/data")
sys.path.append("./Databases/tools")
sys.path.append("./LanguageModel")

import argparse
from Databases import evaluation
from Databases.tools import tools_sequence as sequ
import preprocessing_pca as pca
import generation

import torch
import numpy as np

#from 3DOC.tools import tools_blast as blast #If blasting after sequence generation.
#from 3DOC import output_processing #If running 3DOC pipeline after sequence generation.

if __name__ == "__main__":
    
    #Parses the user input for the arguments --rna_motif and --annotations.
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna_motif", help="RNA-motif, for which an Amino Acid sequence should be generated. Max. 16 entries, recommended: 4 to 8", type=str)
    parser.add_argument('--annotations', help="List of Annotations for Amino Acid Sequence. Please see Databases/annotations_vector_description.csv and go_terms_in_neural_network.txt to see all ontologies. Please separate the onotologies with ', '. Example: 'TAX:Homo, GO:0044822'. The input style for all possible ontologies is: GO-Terms: 'GO:000156'; Taxonomy: 'TAX:Homo', 'TOTALCHARGE-1' (Must be between - 291 and + 200.); Amino-Acid-Types: '#A-004', '#R-004', '#N-004', '#D-004', '#C-004', '#Q-004', '#E-004', '#G-004', '#H-004', '#I-004', '#L-004', '#K-004', '#M-004', '#F-004', '#P-004', '#O-004', '#S-004', '#T-004', '#U-004', '#W-004', '#Y-004', '#V-004', '#B-004', '#Z-004', '#X-004', '#J-004', '#POSITIVE-002', '#NEGATIVE-002', '#POLAR-002', '#NONPOLAR-002', '#HYDROXY-002', '#SULF-002', '#AROMATIC-002'. You must enter an annotation for every category the model has been trained with. ALternatively evaluation.py bool can be adjusted.", type=str)

    args = parser.parse_args()
    rna_motif = args.rna_motif
    annotations = args.annotations
    print(annotations)
    annotations = annotations.split(", ")
    print(annotations)

    #Prepare user-entered RNA-motif, so it can be included in the annotations_vecot.
    rna_motif = sequ.translate_rna_motif_in_anno_vec(rna_motif)
    
    annotations_vector_list = []
    #Create an annotations_vector out of the user-entered annotations and RNA-motif.
    
    annotations_vector, annotations_vector_description, breaking = evaluation.translate_description_to_annotations_vector(annotations, rna_motif)
    np.save("./user_annotations_vector.npy", annotations_vector)    
    if breaking != True:
        annotations_vector = annotations_vector.reshape(1,-1)
        annotations_vector, _ = pca.preprocessing_pca(annotations_vector, annotations_vector)
    np.save("./user_annotations_vector_pca.npy", annotations_vector)

    sequence = generation.generation(annotations_vector)    
    
    fasta_sequence = open("./LanguageModel/output/0.fasta","w")
    fasta_sequence.write(sequence[0])
    fasta_sequence.close()
    
    # If intended, functions can blust against RCSB respectively Swissprot
    #print("Blast against RCSB database.")
    #_, _, _, _, _ = blast.blast_aa_sequence(sequence, "pdb")
    #print("Blast against Swissprot database.")
    #__, _, _, _, _ =blast.blast_aa_sequence(sequence, "swissprot")
    
    # If intended, function can call threedoc_pipline from 3DOC. (Requires running setup.py from tool '3DOC')
    #print("Build PDB with 3DOC and trRosetta.")
    #output_processing.threedoc_pipeline("./LanguageModel/output/", "trRosetta", 3, "best_energy_model")

