from pathlib import Path
from Bio.PDB import *
import numpy as np
import sys
import os
import pyrosetta
import argparse
import subprocess

sys.path.append("./files")
sys.path.append("./tools")
sys.path.append("./trRosetta")

from tools import tools_pyrosetta as toolpyr
from tools import tools_blast as blast

if not os.path.isdir("./output"):
    os.mkdir("./output")

if not os.path.isdir("./output/pdb_files"):
    os.mkdir("./output/pdb_files")
    
if not os.path.isdir("./output/trrosetta"):
    os.mkdir("./output/trrosetta")


def logo():
    print('*********************************************************************')
    print('\
*                ____  _____   ____   _____                         *\n\
*               |___ \|  __ \ / __ \ / ____|                        *\n\
*                 __) | |  | | |  | | |                             *\n\
*                |__ <| |  | | |  | | |                             *\n\
*                ___) | |__| | |__| | |____                         *\n\
*               |____/|_____/ \____/ \_____|                        *')
    print('*                                                                   *')
    print("* Author: Arnoldt, Lucas                                            *")
    print("* iGEM Team Heidelberg 2020                                         *")
    print("* Please email your comments to: igemhd@lucas-arnoldt.de            *")
    print('*********************************************************************')

def threedoc_pipeline(path, program, trrosetta_recurrences, trrosetta_mode):
    '''
    Function enables the execution of the 3DOC pipeline.
    
    :param path:
    :param program:
    :param trrosetta_recurrences:
    :param trrosetta_mode:
    '''
    
    pathlist_fasta = []
    pathlist_seq = []
    breaking = False
    program_saved = program

    for element in Path(path).rglob("*.fasta"):
        pathlist_fasta.append(element)

    pathlist_fasta = np.sort(pathlist_fasta)[::1]

    if len(pathlist_fasta) == 0:
        print("You did not input any FASTA files in the provided directory. Process aborted.")
        breaking = True
    elif program != '3DOC+trRosetta' and program != 'trRosetta':
        print("You did not choose '3DOC+trRosetta' or 'trRosetta' as program for output processing. Please see 'output_processing.py -h'. Process aborted.")
        breaking = True
    elif trrosetta_recurrences > 5:
        print("You did choose more than five recurrences for trRosetta Processing. Please see 'output_processing.py -h'. Process aborted.")
        breaking = True
    elif trrosetta_mode != 'best_energy_model' and trrosetta_mode != 'user_choice':
        print("You did not choose 'best_energy_model' or 'all_ensembles' as mode for trRosetta. Please see 'output_processing.py -h'. Process aborted.")
        breaking = True

    ########################################################################################################################

    if not breaking:
        pyrosetta.init()

        # SAVE LINKERS IN SEPARATE LIST
        protein_dict_fasta = {}
        linker_dict_fasta = {}
        dict_fasta = {}
        sequence_length = 0
        doc_list = []
        trrosetta_list = []
        trrosetta_list_pose_max_energy = []
        trrosetta_list_all_poses = []
        trrosetta_list_all_poses_scores = []
        trrosetta_models_chosen = []

        for sequence_file in pathlist_fasta:
            with open(sequence_file) as fasta_file:

                path = sequence_file.parent
                path = path.joinpath(sequence_file.stem + ".seq")

                with open(path, "w") as seq_file:
                    for line in fasta_file:
                        seq_file.write(line)

                fasta_file.close()
                seq_file.close()

                sequence_length += len(line)

            if sequence_file.stem[1:2] == "-":
                linker_dict_fasta[sequence_file.stem[0:3]] = line
                dict_fasta[sequence_file.stem[0:3]] = line
            else:
                protein_dict_fasta[sequence_file.stem[0:1]] = line
                dict_fasta[sequence_file.stem[0:1]] = line
        
        protein_sequences = list(protein_dict_fasta.values())
        for element in range(0, len(protein_sequences)):
            protein_sequences[element] = protein_sequences[element][:-1]
            
        toolpyr.check_pumby_ppr_length(protein_sequences)

        for key in dict_fasta.keys():
            program = program_saved

            if key[1:2] == "-":
                program = "trRosetta"

            if program == "3DOC+trRosetta":
                doc_list.append(key)
                print("Blasting the amino acid sequence against RCSB database. Blasting may take some time.")
                id, hit, mutated_amino_acids, additional_amino_acids_begin, additional_amino_acids_end = blast.blast_aa_sequence(protein_dict_fasta[key][:-1], "pdb")
                if id != 0:
                    pose = pyrosetta.toolbox.rcsb.pose_from_rcsb(id)
                    blast_fragment_start, blast_fragment_end = blast.find_fragment_start_end(pose, hit)
                    pose = toolpyr.modification_pose(pose, blast_fragment_start, blast_fragment_end, mutated_amino_acids, additional_amino_acids_begin, additional_amino_acids_end)
                    pose.dump_pdb("./output/pdb_files/" + key + ".pdb")
                else:
                    program = "trRosetta"

            if program == "trRosetta":
                trrosetta_list.append(key)

                print("trRosetta Processing started. This may take some time. The usage of a Cluster is recommended.")
                subprocess.call([f"bash trrosetta_pipeline.sh {pathlist_fasta[list(dict_fasta).index(key)].parent} {pathlist_fasta[list(dict_fasta).index(key)].stem} {trrosetta_recurrences}"], shell=True)
                # bash trrosetta_pipeline.sh "/beegfs/home/hd/hd_hd/hd_vu199/trRosetta/" "3_sdgeo"
                print("trRosetta Processing finished.")
         
        pathlist_model = []
        model_scores = []
        
        for element in list(Path("./output/trrosetta").rglob("*.pdb")):
            pathlist_model.append(str(element))
        pathlist_model = np.sort(pathlist_model)[::1]

        trrosetta_list_all_poses = list(pathlist_model)
        counter = 0

        #for protein in range(len(trrosetta_list)):
        if len(pathlist_model) > 0:
            for element in pathlist_model:
                counter += 1
                pose = pyrosetta.io.pose_from_pdb(element)
                score = toolpyr.check_pose_energy(pose)
                model_scores.append(score)
                if counter % trrosetta_recurrences == 0:
                    path_pose_max_energy = pathlist_model[model_scores[counter-trrosetta_recurrences:counter].index(max(model_scores[counter-trrosetta_recurrences:counter]))]
                    trrosetta_list_pose_max_energy.append(path_pose_max_energy)
            trrosetta_list_all_poses_scores.append(model_scores)
            
        for protein in range(len(trrosetta_list)):
            if trrosetta_mode == "best_energy_model":
                trrosetta_models_chosen = trrosetta_list_pose_max_energy
            if trrosetta_mode == "user_choice":
                print(f'For the protein/ linker at the position {trrosetta_list[protein]} the following PDBs with their corresponding energy have been generated with trRosetta (Reason: When starting up output_processing.py trrosetta_recurrences has been set to {trrosetta_recurrences}) :')
                for model in range(0, trrosetta_recurrences):
                    print(f'Protein {trrosetta_list[protein]} - trRosetta model {model}: {trrosetta_list_all_poses[protein*2+model]} - Energy: {trrosetta_list_all_poses_scores[0][protein*2+model]}')
                model_path = input("You have chosen as trrosetta_mode 'user_choice'. Please input the exact path of the model you want to choose to be concatenated in the PDB and input it here. Please dont put the path in quotes.: ")
                trrosetta_models_chosen.append(model_path)

        counter_trrosetta = 0
        
        if len(dict_fasta) > 1:
        
            pose = pyrosetta.io.pose_from_sequence("AAAAAAAA")

            insertion_index = 1
                      
            for element in dict_fasta.keys():
                if element in trrosetta_list:
                    pdb_path = trrosetta_models_chosen[counter_trrosetta]
                    counter_trrosetta += 1
                else:
                    pdb_path = "./output/pdb_files/" + element + ".pdb"
                pose_for_insert = pyrosetta.io.pose_from_pdb(pdb_path)

                mover = pyrosetta.rosetta.protocols.grafting.AnchoredGraftMover(insertion_index, insertion_index+1, pose_for_insert, 0, 0, False)
                mover.apply(pose)
                pose.dump_pdb("./output/3DOC_" + str(insertion_index) + "TEST.pdb")
                insertion_index += len(pose_for_insert.sequence())

            # Delete overhangs, if needed
            #pyrosetta.rosetta.protocols.grafting.delete_region(pose2, len(pose.sequence())-3, len(pose.sequence()))
            #pyrosetta.rosetta.protocols.grafting.delete_region(pose3, 1, 4)

            # Fast Relax
            fr = pyrosetta.rosetta.protocols.relax.FastRelax()
            scorefxn = pyrosetta.rosetta.protocols.loops.get_fa_scorefxn()
            fr.set_scorefxn(scorefxn)
            fr.max_iter(1000)
            if not os.getenv("DEBUG"):
                fr.apply(pose)

            pose.dump_pdb("./output/3DOC_Concatenated.pdb")

            concatenated_sequences = ""

            for key in protein_dict_fasta.keys():
                concatenated_sequences = concatenated_sequences + "_" + key

            aa_sequence = open("./output/3DOC_Conc" + concatenated_sequences + "-PDB.txt","w+")
            aa_sequence.write(pose.sequence())
            aa_sequence.close()

            print("The protein sequences are concatenated and exported as FASTA sequence and model in the output folder.")

        else:

            for element in dict_fasta.keys():
                if element in trrosetta_dict.keys():
                    pdb_path = trrosetta_models_chosen[0]
                else:
                    pdb_path = "./output_tools/pdb_files/" + element + ".pdb"
                pose = pyrosetta.io.pose_from_pdb(pdb_path)

                pose.dump_pdb("./output/3DOC_Concatenated.pdb")

            print("You only inputted one protein sequence. The outputted model can be found in the output folder.")
                      
if __name__ == "__main__":
    logo()
                      
    # BEFEHL: python output_processing.py --path ./files --program geneious_output --trRosetta_recurrences 5 --trRosetta_mode all_ensembles --rna_motif_dict (1,2,3,4,5,6) = AAGGCC

    parser = argparse.ArgumentParser()
    parser.add_argument("--path", help="directory, which contains FASTA files for processing", type=str)
    parser.add_argument("--program", help="'3DOC+trRosetta' or 'trRosetta'", type=str)
    parser.add_argument("--trRosetta_recurrences", help="recommended: '3', type=int", type=int)
    parser.add_argument("--trRosetta_mode", help="'best_energy_model' or 'user_choice'", type=str)
    # parser.add_argument('--rna_motif_dict', help= "dict, of rna motifs for proteins, e.g. {'2,3,4': 'UGAC', '4': 'AGCUCA'}", action = type('', (argparse.Action, ), dict(__call__ = lambda a, p, n, v, o: getattr(n, a.dest).update(dict([v.split('=')])))), default = {}) # anonymously subclassing argparse.Action

    args = parser.parse_args()
    path = args.path
    program = args.program
    trrosetta_recurrences = args.trRosetta_recurrences
    trrosetta_mode = args.trRosetta_mode
    # rna_motif_dict = args.rna_motif_dict
                      
    threedoc_pipeline(path, program, trrosetta_recurrences, trrosetta_mode)

    
