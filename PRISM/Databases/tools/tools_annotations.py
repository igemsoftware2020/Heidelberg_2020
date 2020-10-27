from goatools.semantic import semantic_similarity
import numpy as np

def get_aa_types(sequence):
    """
    Function takes an amino acid sequence as input. It calculates and outputs the number of types of amino acids (e.g. polar, aromatic, ...) and the  total charge of the amino acid sequence.
    Source for amino acid types: Alberts et al. (2015). Molecular Biology of the Cell. 6th Edition. Garland Science, New York.

    :param sequence: amino acid sequence as string
    :return aa_types: number of types of amino acids in sequence
    :return aa_types_description: description of entries in aa_types
    :return total_charge: total charge of amino acid sequence
    """
    
    #Initiate variables
    n_positive = 0
    n_negative = 0
    n_polar = 0
    n_hydroxy = 0
    n_sulf = 0
    n_aromatic = 0
    n_amide = 0

    #Iterate through aminoacids
    for aa in sequence:
        if aa in ['D', 'E']: #negative charge/acidic
            n_negative += 1
            n_polar += 1
        elif aa in ['K', 'H', 'R']: #positive charge/basic
            n_positive += 1
            n_polar += 1
        elif aa in ['F', 'W', 'Y']: #aromatic
            n_aromatic += 1
        elif aa in ['S', 'T', 'Y']: #hydroxy group
            n_hydroxy += 1
            n_polar += 1
        elif aa in ['C', 'M']: #sulfur containing
            n_sulf += 1
        elif aa in ['N', 'Q']: #amide
            n_amide += 1
            n_polar += 1

    n_nonpolar = len(sequence) - n_polar #nonpolar

    total_charge = n_positive - n_negative

    aa_types = [n_positive, n_negative, n_polar, n_nonpolar, n_hydroxy, n_sulf, n_aromatic, n_amide]
    aa_types_description = ['#POSITIVE', '#NEGATIVE', '#POLAR', '#NONPOLAR', '#HYDROXY', '#SULF', '#AROMATIC', '#AMIDE']

    return aa_types, aa_types_description, total_charge

def get_aa_distribution(sequence):
    """
    Function takes an amino acid sequence as input. It calculates and outputs the prevalence of amino acids in the given sequence.

    :param sequence: amino acid sequence as string
    :return distribution: The list distribution includes the prevalence of each amino acid.
    :return distribution_code: The string distribution_code includes all given amino acids.
    """

    distribution_code = 'ARNDCQEGHILKMFPSTWYV'

    distribution = [0 for i in range(0,len(distribution_code))]

    for aa in sequence:
        distribution[distribution_code.find(aa)] += 1

    return distribution, distribution_code

#protein_family_list = ["hello", "sdldg", "sdkg", "hello", "sdldg"]
#param_list = ["0", "1", "2", "3", "4"]

def sorting_entries(protein_family_list, param_list):
    '''
    Function sorts the param_list according to the protein families by using protein_family_list as input.

    :param protein_family_list: list with protein families, which are connected to the entries in param_list
    :param param_list: np array with unsorted sequence matrices or annotations vectors
    :return param_list: np array with sequence matrices or annotations vector sorted by protein families
    '''

    unsorted_list = list(param_list)

    for i in range(0, len(protein_family_list)):
        param_list[i] = unsorted_list[protein_family_list[i]]

    return param_list

def generation_sequence_matrice(sequence, rows):
    """
    Function generates a multi-hot encoding matrice for an amino acid sequence, where a "1" marks the occurence of a specific amino acid.

    :param sequence: amino acid sequence
    :param rows: all possible amino acids in amino acid sequence
    :return sequence_matrice: multi-hot encoding matrice for amino acid sequence
    """

    cols = [0.0] * len(sequence)

    sequence_matrice = np.array([cols] * len(rows))

    for entry in range(0, len(sequence)):
        if sequence[entry] == "0":
            sequence_matrice[rows.index(sequence[entry])][entry] = -float('Inf')
        else: 
            sequence_matrice[rows.index(sequence[entry])][entry] = 1.0

    return sequence_matrice

def translate_sequence_matrice_to_sequence(sequence_matrice, rows):
    """
    Function translates sequence matrice back to a amino acid sequence.

    :param sequence_matrice: multi-hot encoding matrice for amino acid sequence
    :param rows: all possible amino acids in amino acid sequence
    :return sequence: amino acid sequence
    """

    sequence = ""
    
    for element in range(0, len(sequence_matrice[0])):
        for row in range(0, len(rows)):
            if sequence_matrice[row][element] == 1:
                sequence = sequence + rows[row]

    return sequence

def taxonomy_choosing(taxonomy_list, cut_off):
    """
    Function creates list chosen_taxonomies with all taxonomy terms from taxonomy_list, which appear for more than the cut_off number proteins in the record_dict

    :param taxonomy_list: list of all taxonomy terms with their prevalence and frequences
    :param cut_off: number of proteins, which a taxonomy term needs to be appear to be chosen
    :return chosen_taxonomies: list of taxonomies, which appear for more than the cut_off number proteins in the record_dict
    """

    chosen_taxonomies = []

    for taxonomy in taxonomy_list:
        if taxonomy[1][0] > cut_off:
            chosen_taxonomies.append(taxonomy[0])
        else:
            break

    return chosen_taxonomies

def go_terms_choosing(goobj, go_terms_list, go_terms_rbp_list, cut_off):
    '''
    Function is used to choose the GO-Terms used in the annotations_vector by filtering by a cut_off value, the namespace and semantic simialrity.
    
    :param goobj: object with all GO-Terms in Go-Basic file from GO Consortium.
    :param go_terms_list: list of GO-Terms and their prevalence in swissprot database
    :param go_terms_rbp_list: list of GO-Terms, which appear for RBP sequences, and their prevalence in swissprot database
    :param cut_off: number of proteins, which a GO-Term needs to appear for, to not be deleted.
    :return chosen_taxonomies: GO-Terms chosen after applying the cut_off, namespace and semantic similarity dict
    :return go_terms_semantic_similarity_dict: dict with GO-Terms, which are mapped to another GO-Term because of semantic similarity
    '''

    cut_off_go_terms = []
    namespace_go_terms = []
    cut_off_rbp_go_terms = []
    cut_off_rbp = cut_off/len(go_terms_list)*len(go_terms_rbp_list)

    for go_term in go_terms_list:
        if go_term[1][0] > cut_off:
            cut_off_go_terms.append(go_term[0])
        else:
            break

    for go_term in go_terms_rbp_list:
        if go_term[1][0] > cut_off_rbp:
            cut_off_rbp_go_terms.append(go_term[0])
        else:
            break

    for element in cut_off_rbp_go_terms:
        if element not in cut_off_go_terms:
            cut_off_go_terms.append(element)

    for go_term in cut_off_go_terms:
        try:
            if goobj[go_term].namespace == "biological_process" or goobj[go_term].namespace == "molecular_function":
                namespace_go_terms.append(go_term)
        except:
            continue

    level_semantic_similarity_go_terms = []
    level_go_terms = []

    for go_term in namespace_go_terms:
        if goobj[go_term].level < 6:
            level_semantic_similarity_go_terms.append(go_term)

    for go_term in namespace_go_terms:
        if goobj[go_term].level >= 6:
            level_go_terms.append(go_term)

    semantic_similarity_go_terms, go_terms_semantic_similarity_dict = go_terms_semantic_similarity(goobj, level_semantic_similarity_go_terms)

    chosen_taxonomies = level_go_terms + semantic_similarity_go_terms

    return chosen_taxonomies, go_terms_semantic_similarity_dict

def go_terms_children(goobj, go_terms_children_list):
    """
    Function finds and outputs all progeny (children) of a GO-Term list.

    :param goobj: object with all GO-Terms in Go-Basic file from GO Consortium.
    :param go_terms_children_list: list of all GO-Terms, for which all progeny are searched for
    :return go_terms_children_list: list of all GO-Terms with progeny of inputted GO-Terms
    """
    for element in go_terms_children_list:
        try:
            for child in goobj[element].children:
                go_terms_children_list.append(child.id)
        except:
            continue

    return go_terms_children_list

def go_terms_parents(goobj, go_terms_parents_list):
    """
    Function finds and outputs all ancestors (parents) of a GO-Term list.

    :param goobj: object with all GO-Terms in Go-Basic file from GO Consortium.
    :param go_terms_parents_list: list of all GO-Terms, for which all ancestors are searched for
    :return go_terms_parents_list: list of all GO-Terms with ancestors of inputted GO-Terms
    """
    for element in go_terms_parents_list:
        try:
            for parent in goobj[element].parents:
                if parent.id not in go_terms_parents_list:
                    go_terms_parents_list.append(parent.id)
        except:
            continue

    return go_terms_parents_list

def go_terms_semantic_similarity(goobj, go_terms_list):
    """
    Function determines semantic similarity between two GO-terms and assigns via go_terms_semantic_similarity_dict to another go_term.

    :param goobj: object with all GO-Terms in Go-Basic file from GO Consortium.
    :param go_terms_list: list of GO-Terms with a depth below 6.
    :return go_terms_list: list of GO-Terms to which other GO-Terms are assigned for.
    :return go_terms_semantic_similarity_dict: dictionary with assigned GO-Terms to other GO-Terms due to semantic similarity.
    """

    go_terms_semantic_similarity_dict = {}
    semantic_similarity_go_terms_list = []

    for i in range(0, len(go_terms_list)-1):
        for j in range(i+1, len(go_terms_list)):
            if go_terms_list[j] not in semantic_similarity_go_terms_list and go_terms_list[i] not in semantic_similarity_go_terms_list:
                similarity = semantic_similarity(go_terms_list[i], go_terms_list[j], goobj)
                if type(similarity) is float:
                    if similarity > 0.8:
                        go_terms_semantic_similarity_dict[go_terms_list[j]] = go_terms_list[i]
                        semantic_similarity_go_terms_list.append(go_terms_list[j])

    for i in range(len(semantic_similarity_go_terms_list)-1, -1, -1):
        go_terms_list.pop(go_terms_list.index(semantic_similarity_go_terms_list[i]))

    return go_terms_list, go_terms_semantic_similarity_dict


if __name__ == '__main__':
    s = input()
    print(get_aa_distribution (s))
