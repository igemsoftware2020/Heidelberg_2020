import rna_denovo_preparation as rna
import pyrosetta
import argparse
import os
    
def protein_rna_sequence(sequence_protein, sequence_rna):
    '''
    Function prepares the file "protein_sequence_concatenated_with_rna.tx" with a combined sequence of protein sequence and RNA-motif for usage in Rosetta RNP Structure modelling.
    
    :param sequence_protein: sequence of protein
    :param sequence_rna: RNA-motif
    '''
    
    sequence_protein_rna = sequence_protein + sequence_rna
    
    file_protein_sequence_concatenated_with_rna = open("./rnp_prediction/protein_sequence_concatenated_with_rna.txt", "w") 
    file_protein_sequence_concatenated_with_rna.write(sequence_protein_rna)
    file_protein_sequence_concatenated_with_rna.close() 
    
def protein_rna_secstruct(secstruct_rna, sequence_protein, sequence_rna):
    '''
    Function prepares the file "protein_rna_secstruct.txt" with a combined secondary structure of the RNA-motif and the protein for usage in Rosetta RNP Structure modelling.
    
    :param secstruct_rna: Secondary structure of RNA-motif
    :param sequence_protein: sequence of protein
    :param sequence_rna: sequence of RNA
    '''
    
    secstruct_protein = len(sequence_protein) * "."
    
    protein_rna_secstruct = secstruct_protein + secstruct_protein
    sequence_protein_rna = sequence_protei + sequence_rna
    
    file_protein_rna_secstruct = open("./rnp_prediction/protein_rna_secstruct.txt","w") 
    file_protein_rna_secstruct.write(protein_rna_secstruct)
    file_protein_rna_secstruct.write(sequence_protein_rna)
    file_protein_rna_secstruct.close() 
    
def protein_rna_pdb(pdb_protein, pdb_rna):
    '''
    Function combines a Protein and a RNA PDB to the single pdb file "unbound_protein_and_RNA.pdb"for usage in Rosetta RNP Structure modelling.
    
    :param pdb_protein: path to the PDB file of the protein which should be combined with the PDB of the RNA
    :param pdb_rna: path to the PDB file of the RNA which should be combined with the PDB of the Protein
    '''
    
    file_pdb_protein = open(pdb_protein, "r")
    file_pdb_rna = open(pdb_rna, "r")
   
    file_pdb_protein_readlines = file_pdb_protein.readlines()
    file_pdb_rna_readlines = file_pdb_rna.readlines()
    
    #NEW FILE
    file_unbound_protein_and_RNA = open("./rnp_prediction/unbound_protein_and_RNA.pdb", "w")
    
    
    for line_pdb in file_pdb_protein_readlines:
        if line_pdb[0:4] == "ATOM":
            file_unbound_protein_and_RNA.write(line_pdb)
            line_pdb_split = line_pdb.split(" ")
   
    line_pdb_split_processed = []
    for element in line_pdb_split:
        if element != "":
            line_pdb_split_processed.append(element)
            
    position_number_pdb = int(line_pdb_split_processed[1])
    base_number_pdb = int(line_pdb_split_processed[5])
    
    count_rna = 0
    base_number_rna_entry_before = -100
    count_base_number = 0
    position_before = position_number_pdb
    
    for line_rna in file_pdb_rna_readlines:
        count_rna += 1
        if line_rna[0:6] == "HETATM":
            position = str(position_before + 1)
            position_string = " " * (7-len(position)) + position
            base_number_rna = int(line_rna[22:26])
            if base_number_rna != base_number_rna_entry_before: 
                base_number = str(base_number_pdb + 3 + count_base_number)
                count_base_number += 1
                base_number_string = " " * (4-len(base_number)) + base_number
            base_number_rna_entry_before = base_number_rna
            line_rna_processed = "ATOM" + position_string + line_rna[11:22] + base_number_string + line_rna[26:]
            file_unbound_protein_and_RNA.write(line_rna_processed)
            position_before = int(position)

    file_unbound_protein_and_RNA.write("END")
    
    file_pdb_protein.close()
    file_pdb_rna.close()
    file_unbound_protein_and_RNA.close()   
	
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--secstruct_rna", help = "Secondary structure of RNA modelled in RNP-Complex. Optional. Otherwise calculated.", type = str) #OPTIONAL
    parser.add_argument("--pdb_protein", help = "PDB of Protein, which is modelled in RNP-Complex.", type = str, required=True)
    parser.add_argument("--pdb_rna", help = "PDB of RNA, which is modelled in RNP-Complex.", type = str, required=True)
	
    args = parser.parse_args()
    secstruct_rna = args.secstruct_rna
    pdb_protein = args.pdb_protein
    pdb_rna = args.pdb_rna
    
    if not os.path.isdir("rnp_prediction"):
        os.mkdir("rnp_prediction")
    
    pyrosetta.init()
    
    pose_protein = pyrosetta.io.pose_from_pdb(pdb_protein)
    pose_rna = pyrosetta.io.pose_from_pdb(pdb_rna)

    if not (args.secstruct_rna):
        rna_motif = pose_rna.sequence()
        secstruct_rna = rna.secondary_structure_rna(rna_motif)
        print("Secondary structure RNA:")
        print(secstruct_rna)
    
    sequence_protein = pose_protein.sequence()
    sequence_rna = pose_rna.sequence()

    protein_rna_sequence(sequence_protein, sequence_rna)
    protein_rna_secstruct(secstruct_rna, sequence_protein, sequence_rna)
    protein_rna_pdb(pdb_protein, pdb_rna)
