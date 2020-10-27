from Bio import SeqIO
from tools import tools_annotations as anno
from tools import tools_sequence as sequ

def extraction_go_term_taxonomy(record_dict, record):
    """
    Function extracts the go terms and taxonomy of a specific record in record_dict.

    :param record_dict: dictionary of all records of swissprot
    :param record: uniprot accession number of record in record dict
    :return go_terms: go terms of record
    :return taxonomy: taxonomy of record
    """
    record = record_dict[record]

    taxonomy = record.annotations['taxonomy'] #LISTE
    go_terms = extraction_go_terms(record.dbxrefs)

    return go_terms, taxonomy

def extraction_annotations(record_dict, record, unified_record_sequence_length):
    """
    Function extracts several annotations of a specific record in record_dict.

    :param record_dict: dictionary of all records of swissprot
    :param record: uniprot accession number of record in record dict
    :param unified_record_sequence_length: the unified length of all sequences used in the neural network
    :return id: uniprot accession number of record
    :return sequence_list: sequences of record
    :return mask_list: masks for sequences of record
    :return taxonomy: taxonomy of record
    :return ncbi_taxid: ncbi_taxid of record
    :return go_terms: go terms of record
    :return family: family of record
    """
    record = record_dict[record]

    id = record.id
    sequence = str(record.seq)
    sequence = sequ.resolve_amino_acids(sequence)
    sequence_list, mask_list = sequ.cut_sequence_length(sequence, unified_record_sequence_length)
    taxonomy = record.annotations['taxonomy']
    ncbi_taxid = record.annotations['ncbi_taxid'][0]
    go_terms = extraction_go_terms(record.dbxrefs)
    family = extraction_family(record.dbxrefs)

    return id, sequence_list, mask_list, taxonomy, ncbi_taxid, go_terms, family

def extraction_family(dbxrefs):
    """
    Function extracts pfam family from description dbxrefs for a specific sequence.

    :param dbxrefs: string containing annotations for a sequence such as go terms
    :return family: pfam family in description dbxrefs for a specific sequence
    """
    family_record = ""

    for element in range(0, len(dbxrefs)):
        if dbxrefs[element][0:5] == "Pfam:":
            family_record = dbxrefs[element][5:]
            break

    return family_record

def extraction_go_terms(dbxrefs):
    """
    Function extracts GO-Terms from description dbxrefs for a specific sequence.

    :param dbxrefs: string containing annotations for a sequence such as go terms
    :return go_terms: list of go terms in description dbxrefs for a specific sequence
    """
    go_terms = []

    for element in range(0, len(dbxrefs)):
        if dbxrefs[element][0:3] == "GO:":
            go_terms.append(dbxrefs[element][3:])

    return go_terms


def annotations_frequency(terms_dict, counter):
    """
    Function calculates frequency of GO-Terms in go_terms_dict and outputs the go_terms_list.

    :param terms_dict: dictionary with go terms and their prevalence
    :param counter: number of GO-Terms in total
    :return terms_list: list with GO-Terms and their prevalence and frequency
    """
    family = 0

    # Calculation go term frequency in percent
    for term in terms_dict:
        freq = terms_dict[term][0] / counter
        terms_dict[term].append(freq)

    terms_list = list(terms_dict.items())  # "Conserving" the dictionary
    terms_list.sort(key=lambda x: x[1][1], reverse=True)  # Sorting by frequency

    return terms_list
