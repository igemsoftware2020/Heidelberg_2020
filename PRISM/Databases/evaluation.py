import os
import sys
sys.path.append(os.getcwd())
sys.path.append("./data")
sys.path.append("./tools")
sys.path.append("./Databases")
sys.path.append("./Databases/data")
sys.path.append("./Databases/tools")
sys.path.append("./LanguageModel")

from Databases.tools import tools_annotation_vector as anno_vec

import pickle
import numpy as np

def translate_description_to_annotations_vector(user_input_annotations, rna_motif):
    """
    Function can be called with the annotations provided by the user for creation of a new RBP sequence by the neural network. The function then creates the annotations_vector and annotations_vector_description for the user input.

    :param user_input_annotations: user input for annotations for creation a new RBP sequence out of the network
    :return annotations_vector: annotations vector with the rows representing the annotations provided by the user
    :return annotations_vector_description: vector with descriptions for the rows in the annotations_vector
    :return breaking: bool, whether an annotation was added, which does not exist
    """

    go_terms_annotation = False
    taxonomy_annotation = False
    total_charge_annotation = False
    aa_distribution_annotation = False
    aa_types_annotation = False
    ligand_binding_annotation = False

    # INPUT GOTERM: e.g. "GO:000156"
    # INPUT TAXONOMY: TBD, maybe "TAX:Homo"
    # INPUT TOTAL CHARGE TYPE: e.g. "TOTALCHARGE-1"      #DETERMINE CHARGE RANGE
    # INPUT AA DISTRIBUTION CHARGE: e.g. "POSITIVE-003-NEGATIVE-234"
    # INPUT AA DISTRIBUTION POLARITY: e.g. "POLAR-003-UNPOLAR-234"
    # INPUT AA DISTRIBUTION POLARITY: e.g. "LIGAND-Y"/ "LIGAND-N"

    go_terms = []
    taxonomy = []
    total_charge = 0
    distribution_code = 'ARNDCQEGHILKMFPOSTWYVBZXJ'
    aa_distribution = [0 for i in range(0,len(distribution_code))]
    aa_types = [0, 0, 0, 0, 0, 0, 0, 0]
    aa_types_description = ["#POSITIVE", "#NEGATIVE", "#POLAR", "#NONPOLAR", "#HYDROXY", "#SULF", "#AROMATIC", "#AMIDE"]
    id = 0
    breaking = False

    go_terms_semantic_similarity_dict_file = open("./Databases/data/swissprot/go_terms_semantic_similarity_dict.pkl", "rb")
    go_terms_semantic_similarity_dict = pickle.load(go_terms_semantic_similarity_dict_file)
    go_terms_semantic_similarity_dict_file.close()
    bindingdb_uniprot_ids = list(np.load("./Databases/data/binding_db/bindingdb_uniprot_ids.npy"))
    go_terms_selection = list(np.load("./Databases/data/swissprot/go_terms_annotations_vector_pre_compiled.npy"))
    taxonomy_selection = list(np.load("./Databases/data/swissprot/taxonomy_annotations_vector_pre_compiled.npy"))

    for element in range(0, len(user_input_annotations)):
        if user_input_annotations[element][0:2] == "GO":
            go_terms_annotation = True
            go_terms.append(user_input_annotations[element])
        elif user_input_annotations[element][0:3] == "TAX":
            taxonomy_annotation = True
            taxonomy.append(user_input_annotations[element])
        elif user_input_annotations[element][0:11] == "TOTALCHARGE":
            total_charge_annotation = True
            total_charge = float(user_input_annotations[element][13:])
        elif user_input_annotations[element][0:1] == "#" and user_input_annotations[element][2:3] == "-":
            aa_distribution_annotation = True
            prevalence = float(user_input_annotations[element][3:])
            aa_distribution[distribution_code.find(user_input_annotations[element][1:2])] = prevalence
        elif user_input_annotations[element][0:1] == "#":
            aa_types_annotation = True
            aa_type = (user_input_annotations[element].split("-"))[0]
            prevalence = float((user_input_annotations[element].split("-"))[1])
            aa_types[aa_types_description.index(aa_type)] = prevalence
        elif user_input_annotations[element][0:6] == "LIGAND":
            ligand_binding_annotation = True
            id = bindingdb_uniprot_ids[0]
        else:
            print("You have added an annotation, which does not exist in the Neural Network. Please check the given annotations list.")
            breaking = True
            break

    rna_motif_blank = [0.0] * 4 * 12
    annotations_vector, annotations_vector_description = anno_vec.creation_blank_annotations_vector(go_terms_annotation, go_terms_selection, taxonomy_annotation, taxonomy_selection, total_charge_annotation, aa_distribution_annotation, aa_types_annotation, ligand_binding_annotation, rna_motif_blank)
    
    annotations_vector = anno_vec.entry_specific_modification_annotation_vector(id, annotations_vector, annotations_vector_description, go_terms, go_terms_semantic_similarity_dict, taxonomy, total_charge, aa_distribution, distribution_code, aa_types, aa_types_description, bindingdb_uniprot_ids, rna_motif)
   
    return annotations_vector, annotations_vector_description, breaking


if __name__ == '__main__':
    user_input_annotations = input("Please type in the annotations you regard as necessary for a new RBP as an array in the type (STRING1, STRING2, ...). Please spell/ write the annotations exactly like the given annotations table.: ")
    annotations_vector, annotations_vector_description = translate_description_to_annotations_vector(user_input_annotations)
