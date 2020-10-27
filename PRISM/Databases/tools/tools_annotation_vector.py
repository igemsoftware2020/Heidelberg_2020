import numpy as np
from tools import tools_normalisation as norm
from tools import tools_binding_db as db

def creation_blank_annotations_vector(go_terms, go_terms_selection, taxonomy, taxonomy_selection, total_charge, aa_distribution, aa_types, ligand_binding, rna_motif):
    """
    The function creates a blank annotations_vector with a another vector of the same length (annotations_vector_description) with the description for every row of the annotations vector. The annotations_vector contains rows for all chosen annotations.

    :param go_terms: True/ False to determine, whether the blank annotations vector should provide rows for go_terms
    :param go_terms_selection: selection of go terms of all records in swissprot to use in annotations vectors
    :param taxonomy: True/ False to determine, whether the blank annotations vector should provide rows for taxonomy
    :param taxonomy_selection: selection of taxonomy terms of all record in swissprot to use in annotations vectors
    :param total_charge: True/ False to determine, whether the blank annotations vector should provide rows for the total charge of aa sequence
    :param aa_distribution: True/ False to determine, whether the blank annotations vector should provide rows for the distribution of amino acids in sequence
    :param aa_types: True/ False to determine, whether the blank annotations vector should provide rows for the prevalence of types of aa (polar, positive, ...)
    :param ligand_binding: True/ False to determine, whether the blank annotations vector should provide rows for the characteristic of ligand-binding (Yes/ No)
    :param rna_motif: RNA-motif bound by the protein in the record in swissprot
    :return annotations_vector: blank annotations vector with the rows representing the params
    :return annotations_vector_description: vector with descriptions for the rows in the blank annotations_vector
    """

    go_terms_vector = []
    taxonomy_vector = []
    total_charge_vector = []
    aa_distribution_vector = []
    aa_types_vector = []
    ligand_binding_vector = []
    rna_motif_vector = []
    spacer = [0] * 48

    for i in range(0, len(taxonomy_selection)):
        modified_taxonomy = 'TAX:' + taxonomy_selection[i]
        taxonomy_selection[i] = modified_taxonomy
    if go_terms == True:
        go_terms_vector = list(go_terms_selection)
    if taxonomy == True:
        taxonomy_vector = list(taxonomy_selection)
    if total_charge == True:
        total_charge_vector = ["Total Charge[Between -291 and +200]"]
    if aa_distribution == True:
        aa_distribution_vector = ['#A', '#R', '#N', '#D', '#C', '#Q', '#E', '#G', '#H', '#I', '#L', '#K', '#M', '#F', '#P', '#S', '#T', '#W', '#Y', '#V']
    if aa_types == True:
        aa_types_vector = ['#positive', '#negative', '#polar', '#nonpolar', '#hydroxy', '#sulf', '#aromatic', '#amide']
    if ligand_binding == True:
        ligand_binding_vector = ["Ligand Binding [Yes = 1, No = 0]"]
    if rna_motif == True:
        rna_motif_vector = ['ResPos1-A', 'ResPos1-C', 'ResPos1-G', 'ResPos1-U', 'ResPos2-A', 'ResPos2-C', 'ResPos2-G', 'ResPos2-U', 'ResPos3-A', 'ResPos3-C', 'ResPos3-G', 'ResPos3-U', 'ResPos4-A', 'ResPos4-C', 'ResPos4-G', 'ResPos4-U', 'ResPos5-A', 'ResPos5-C', 'ResPos5-G', 'ResPos5-U', 'ResPos6-A', 'ResPos6-C', 'ResPos6-G', 'ResPos6-U', 'ResPos7-A', 'ResPos7-C', 'ResPos7-G', 'ResPos7-U', 'ResPos8-A', 'ResPos8-C', 'ResPos8-G', 'ResPos8-U', 'ResPos9-A', 'ResPos9-C', 'ResPos9-G', 'ResPos9-U', 'ResPos10-A', 'ResPos10-C', 'ResPos10-G', 'ResPos10-U', 'ResPos11-A', 'ResPos11-C', 'ResPos11-G', 'ResPos11-U', 'ResPos12-A', 'ResPos12-C', 'ResPos12-G', 'ResPos12-U']

    length_annotations_vector = len(go_terms_vector) + len(taxonomy_vector) + len(total_charge_vector) + len(aa_distribution_vector) + len(aa_types_vector) + len(ligand_binding_vector) + len(spacer)
    annotations_vector = np.array([0.0] * length_annotations_vector, dtype=float)
    annotations_vector_description = go_terms_vector + taxonomy_vector + total_charge_vector + aa_distribution_vector + aa_types_vector + ligand_binding_vector + spacer + rna_motif_vector

    try:
        with open('./data/annotations_vector_description.csv', 'w') as filehandle:
            for listitem in annotations_vector_description:
                filehandle.write('%s\n' % listitem)
    except:
        pass
                
    annotations_vector_description_array = np.asarray(annotations_vector_description)
    try:
        np.save("./data/annotations_vector_description.npy", annotations_vector_description_array)
    except:
        pass

    return annotations_vector, annotations_vector_description


def entry_specific_modification_annotation_vector(id, annotations_vector, annotations_vector_description, go_terms, go_terms_semantic_similarity_dict, taxonomy, total_charge, aa_distribution, distribution_code, aa_types, aa_types_description, bindingdb_uniprot_ids, rna_motif):
    """
    Function modifies the annotation vector based on the record information.

    :param id: uniprot accession code
    :param annotations_vector: blank annotations vector with the rows representing the params
    :param annotations_vector_description: vector with descriptions for the rows in the blank annotations_vector
    :param go_terms: go terms of the record
    :param go_terms_semantic_similarity_dict: dictionary of GO-terms which are assigned to each other because of semantic similarity
    :param taxonomy: taxonomy information of the record
    :param total_charge: total charge of the sequence of the record
    :param aa_distribution: prevalence of amino acids in the sequence of the record
    :param distribution_code: description for information given in aa_distribution
    :param aa_types: prevalence of amino acids of specific types (positive, polar, ...)
    :param aa_types_description: description for information given in aa_types
    :param bindingdb_uniprot_ids: information (Y/ N) whether sequence of the record is able to bind a ligand or not
    :param rna_motif: RNA-motif bound by the protein in the record in swissprot
    :return annotations_vector: annotations vector with the rows representing the params inputted in this function
    """

    # Replace
    for i in range(0, len(go_terms)):
        try:
            go_term_semantic_similarity = go_terms_semantic_similarity_dict[go_terms[i]]
            go_terms[i] = go_term_semantic_similarity
        except:
            continue

    for element in range(0, len(annotations_vector_description)):
        if annotations_vector_description[element] != 0:
            if annotations_vector_description[element][0:2] == "GO":
                go_terms_bool = annotations_vector_description[element] in go_terms
                if go_terms_bool == True:
                    annotations_vector[element] = 1.0
            elif annotations_vector_description[element][0:4] == "TAX:":
                taxonomy_bool = annotations_vector_description[element][4:] in taxonomy
                if taxonomy_bool == True:
                    annotations_vector[element] = 1.0
            elif annotations_vector_description[element] == 'Total Charge[Between -291 and +200]':
                annotations_vector[element] = norm.normalisation_radial_basis_function(total_charge, 4.238972298883644, 136.56931758567416)
            elif annotations_vector_description[element][0:1] == "#":
                if annotations_vector_description[element] in aa_types_description:
                    annotations_vector[element] = norm.normalisation_radial_basis_function(aa_types[aa_types_description.index(annotations_vector_description[element])], 519.3488124287404, 76491.78050523378)
                elif annotations_vector_description[element][1:1] in distribution_code:
                    annotations_vector[element] = norm.normalisation_radial_basis_function(aa_distribution[distribution_code.index(annotations_vector_description[element][1:1])], 682.1225807253651, 127911.0850291213)
            elif annotations_vector_description[element] == 'Ligand Binding [Yes = 1, No = 0]':
                if id in bindingdb_uniprot_ids:
                    annotations_vector[element] = 1.0
            else:
                binding_db_bool = db.update_ligand_binding(id, bindingdb_uniprot_ids)
                if binding_db_bool == True:
                    annotations_vector[element] = 1.0
            
    annotations_vector = np.concatenate((annotations_vector, rna_motif))

    return annotations_vector

def translate_annotations_vector_to_description(annotations_vector, annotations_vector_description):
    """
    Function translates an annotations_vector for a specific record to the description.

    :param annotations_vector: annotations vector with the rows representing the annotations provided by the user
    :param annotations_vector_description: vector with descriptions for the rows in the annotations_vector
    """
    for row in range(0, len(annotations_vector)):
        if annotations_vector[row] == 1:
            print(annotations_vector_description[row])
