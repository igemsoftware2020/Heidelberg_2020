import os
import sys
sys.path.append(os.getcwd())
sys.path.append("./data")
sys.path.append("./tools")
sys.path.append("./Databases")
sys.path.append("./Databases/data")

from ftplib import FTP
import gzip
import shutil
import urllib.request
import zipfile
import pickle
from Bio import SeqIO
import numpy as np
import pandas as pd
from goatools.obo_parser import GODag

from tools import tools_annotation_vector as anno_vec
from tools import tools_swissprot as swiss
from tools import tools_binding_db as db
from tools import tools_attract as att
from tools import tools_annotations as anno
from tools import tools_sequence as sequ

#Creation Folders for Data and Automatic Download of Databases, which are used in the Neural Network
print("Folder creation and database download starting.")

if not os.path.isdir("./data"):
    os.mkdir("./data")

if not os.path.isdir("./data/swissprot"):
    os.mkdir("./data/swissprot")

    ftp = FTP('ftp.uniprot.org')
    ftp.login()
    ftp.cwd('pub/databases/uniprot/current_release/knowledgebase/complete')
    ftp.retrlines('LIST')
    with open('uniprot_sprot.dat.gz', 'wb') as fp:
        #ftp.retrbinary('RETR uniprot_sprot.dat.gz', fp.write)
        ftp.retrbinary('RETR uniprot_sprot.dat.gz', open(os.path.join("./data/swissprot", "uniprot_sprot.dat.gz"), "wb").write)
    ftp.quit()

    with gzip.open('./data/swissprot/uniprot_sprot.dat.gz', 'rb') as f_in:
        with open('./data/swissprot/uniprot_sprot.dat', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

if not os.path.isdir("./data/binding_db"):
    os.mkdir("./data/binding_db")

    binding_db_url = 'https://www.bindingdb.org/bind/BindingDB_UniProt.txt'
    urllib.request.urlretrieve(binding_db_url, './data/binding_db/BindingDB_UniProt.txt')

if not os.path.isdir("./data/attract"):
    os.mkdir("./data/attract")

    attract_url = 'https://attract.cnic.es/attract/static/ATtRACT.zip'
    urllib.request.urlretrieve(attract_url, './data/attract/ATtRACT.zip')

    with zipfile.ZipFile('./data/attract/ATtRACT.zip', 'r') as zip_ref:
        zip_ref.extractall('./data/attract')

if not os.path.isdir("./data/go"):
    os.mkdir("./data/go")

    go_url = 'http://current.geneontology.org/ontology/go-basic.obo'
    urllib.request.urlretrieve(go_url, './data/go/go-basic.obo')

if not os.path.isdir("./data/training"):
    os.mkdir("./data/training")

if not os.path.isdir("./data/finetuning"):
    os.mkdir("./data/finetuning")


#######################################################################################################################

# Loading record_dict which contains all records from swissprot database

record_dict = SeqIO.index("./data/swissprot/uniprot_sprot.dat", "swiss")

#######################################################################################################################

# Counting GO-Terms and taxonomies in record_dict and deciding for go_terms and taxonomies to use for annotations for the neural network.
print("Determining GO-Terms and taxonomies for neural network starting.")

# Alternative: Loading pre-calculated taxonomies, GO-terms list and go_terms_semantic_similarity_dict to use for annotations:
# go_terms_selection = np.load("./data/swissprot/go_terms_annotations_vector_pre_compiled.npy")
# taxonomy_selection = np.load("./data/swissprot/taxonomy_annotations_vector_pre_compiled.npy")
# go_terms_selection = list(go_terms_selection)
# taxonomy_selection = list(taxonomy_selection)
# go_terms_semantic_similarity_dict_file = open("./data/swissprot/go_terms_semantic_similarity_dict.pkl", "rb")
# go_terms_semantic_similarity_dict = pickle.load(go_terms_semantic_similarity_dict_file)
# go_terms_semantic_similarity_dict_file.close()

# Loading of go-basic.obo file with list of GO-Terms from GO Consortium
goobj = GODag("./data/go/go-basic.obo", optional_attrs='relationship')

# Definition of GO-Terms, which stand for catalytic activity and finding all progeny (children) of them.
go_terms_catalytic_activity_list = ["GO:0140096", "GO:0140097", "GO:0140098"]
go_terms_catalytic_activity_list = anno.go_terms_children(goobj, go_terms_catalytic_activity_list)

# Extraction of all GO-Terms and Taxonomy from records in record dict for all entries and RBPs for GO-Terms only.

go_terms_dict = {}
go_terms_rbp_dict = {}
taxonomy_dict = {}
go_terms_counter = 0
rbp_go_terms_counter = 0
taxonomy_counter = 0

for record in record_dict:

    go_terms, taxonomies = swiss.extraction_go_term_taxonomy(record_dict, record)

    counter = 0

    for go_term in go_terms:
        if go_term not in go_terms_catalytic_activity_list:
            counter += 1

    if counter == len(go_terms):

        go_terms = anno.go_terms_parents(goobj, go_terms)

        if "GO:0003723" in go_terms or "GO:0044822" in go_terms or "GO:0000498" in go_terms:
            for go_term in go_terms:
                if go_term not in go_terms_rbp_dict:
                    go_terms_rbp_dict[go_term] = [1]
                else:
                    go_terms_rbp_dict[go_term][0] += 1
                rbp_go_terms_counter += 1

        for go_term in go_terms:
            if go_term not in go_terms_dict:
                go_terms_dict[go_term] = [1]
            else:
                go_terms_dict[go_term][0] += 1
            go_terms_counter += 1

        for taxonomy in taxonomies:
            if taxonomy not in taxonomy_dict:
                taxonomy_dict[taxonomy] = [1]
            else:
                taxonomy_dict[taxonomy][0] += 1
            taxonomy_counter += 1

# Calculation of frequency of GO-Terms and Taxonomy in contrast to total terms.

go_terms_list = swiss.annotations_frequency(go_terms_dict, go_terms_counter)
go_terms_rbp_list = swiss.annotations_frequency(go_terms_rbp_dict, rbp_go_terms_counter)
taxonomy_list = swiss.annotations_frequency(taxonomy_dict, taxonomy_counter)

# Function chooses all GO-Terms, which fulfill certain criteria (Appear for more than the cut_off number proteins in the record_dict,
cut_off = 1000
go_terms_selection, go_terms_semantic_similarity_dict = anno.go_terms_choosing(goobj, go_terms_list, go_terms_rbp_list, cut_off)
go_terms_selection = go_terms_selection[0:904]

# Function saves all taxonomies which appear for more than the cut_off number proteins in the record_dict
cut_off = 7000
taxonomy_selection = anno.taxonomy_choosing(taxonomy_list, cut_off)

# Saving go terms and taxonomy, which are used in the annotations vector and the neural network

selected_go_terms = open('./data/go_terms_in_neural_network.txt', 'w')
selected_go_terms.writelines('GO-Terms used in the Neural Network:' + '\n')
for go_term_selected in go_terms_selection:
    selected_go_terms.writelines(go_term_selected + "-" + goobj[go_term_selected].name + '\n')
selected_go_terms.close()

go_terms_selection_array = np.asarray(go_terms_selection)
np.save("./data/swissprot/go_terms_annotations_vector_pre_compiled.npy", go_terms_selection_array)

selected_taxonomy = open('./data/taxonomies_in_neural_network.txt', 'w')
selected_taxonomy.writelines('Taxonomies used in the Neural Network:' + '\n')
for taxonomy_selected in taxonomy_selection:
    selected_taxonomy.writelines(taxonomy_selected + '\n')
selected_taxonomy.close()

taxonomy_selection_array = np.asarray(taxonomy_selection)
np.save("./data/swissprot/taxonomy_annotations_vector_pre_compiled.npy", taxonomy_selection_array)

go_terms_semantic_similarity_dict_file = open("./data/swissprot/go_terms_semantic_similarity_dict.pkl", "wb")
pickle.dump(go_terms_semantic_similarity_dict, go_terms_semantic_similarity_dict_file)
go_terms_semantic_similarity_dict_file.close()

#######################################################################################################################

# Choosing Annotations, which should be used in the Neural Network (go-terms, taxonomy, total_charge,
# aa_distribution, aa_types, ligand_binding)
print("Iterating through Swissprot and creation of annotations vectors starting.")

go_terms_annotation = True
taxonomy_annotation = True
total_charge_annotation = True
aa_distribution_annotation = True
aa_types_annotation = True
ligand_binding_annotation = True
rna_motif_annotation = True

# Creation of a blank annotations vector with the right length to fit all annotations chosen above.
annotations_vector, annotations_vector_description = anno_vec.creation_blank_annotations_vector(go_terms_annotation, go_terms_selection, taxonomy_annotation, taxonomy_selection, total_charge_annotation, aa_distribution_annotation, aa_types_annotation, ligand_binding_annotation, rna_motif_annotation)
annotations_vector_saved = np.array(annotations_vector)

# Extraction of all uniprot ids binding a ligand from the binding db database
bindingdb_uniprot_ids = db.extraction_uniprot_ids_bindingdb()
bindingdb_uniprot_ids = np.asarray(bindingdb_uniprot_ids)
np.save("./data/binding_db/bindingdb_uniprot_ids.npy", bindingdb_uniprot_ids)
# Alternatively load with:
# bindingdb_uniprot_ids = np.load("./data/binding_db/bindingdb_uniprot_ids.npy")

# Iterating through the record_dict containing all records of swissprot to extract the annotations and the sequence
# and to create an annotations vector

annotations_vector_list = []
sequences_list = []
mask_list = []
protein_family_list = []
unified_record_sequence_length = 512
rna_motif = ""
rna_motif = sequ.translate_rna_motif_in_anno_vec(rna_motif)

for record in record_dict:

    annotations_vector = np.array(annotations_vector_saved, dtype=float)
    id, sequence_list, mask_list_sequence_list, taxonomy, ncbi_taxid, go_terms, family = swiss.extraction_annotations(record_dict, record, unified_record_sequence_length)

    counter = 0

    for go_term in go_terms:
        if go_term not in go_terms_catalytic_activity_list:
            counter += 1

    if counter == len(go_terms):

        for sequence, mask in zip(sequence_list, mask_list_sequence_list):

            aa_distribution, distribution_code = anno.get_aa_distribution(sequence)
            aa_types, aa_types_description, total_charge = anno.get_aa_types(sequence)

            # Modification of blank annotations vector based on annotations of the record extracted from the record_dict
            annotations_vector = anno_vec.entry_specific_modification_annotation_vector(id, annotations_vector, annotations_vector_description, go_terms, go_terms_semantic_similarity_dict, taxonomy, total_charge, aa_distribution, distribution_code, aa_types, aa_types_description, bindingdb_uniprot_ids, rna_motif)
            rows_amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '0']
            #sequence = anno.generation_sequence_matrice(sequence, rows_amino_acids)
            sequence = sequ.generation_sequence_one_hot_preparation(sequence, rows_amino_acids)

            if len(sequence) == unified_record_sequence_length:
                if annotations_vector_list is None:
                    annotations_vector_list = annotations_vector[None]
                else:
                    annotations_vector_list.append(annotations_vector)
                
                sequences_list.append(sequence)
                mask_list.append(mask)
                protein_family_list.append(family)
        
            annotations_vector = np.array(annotations_vector_saved, dtype=float)

# Sequences and annotations vectors are sorted for the protein families.
protein_family_list_argsort = np.argsort(protein_family_list)
protein_family_list = anno.sorting_entries(protein_family_list_argsort, protein_family_list)
sequences_list = anno.sorting_entries(protein_family_list_argsort, sequences_list)
mask_list = anno.sorting_entries(protein_family_list_argsort, mask_list)
annotations_vector_list = anno.sorting_entries(protein_family_list_argsort, annotations_vector_list)

protein_family_list = np.asarray(protein_family_list)
sequences_list = np.asarray(sequences_list)
mask_list = np.asarray(mask_list)
annotations_vector_list = np.asarray(annotations_vector_list, dtype=float)

# Saving the rna_motif, binding sequences and scores for use in the neural network
np.save("/beegfs/work/ws/hd_vu199-preprocessing_arnoldt-0/data/training/protein_families.npy", protein_family_list)
np.save("/beegfs/work/ws/hd_vu199-preprocessing_arnoldt-0/data/training/masks.npy", mask_list)
np.save("/beegfs/work/ws/hd_vu199-preprocessing_arnoldt-0/data/training/sequences.npy", sequences_list)
np.save("/beegfs/work/ws/hd_vu199-preprocessing_arnoldt-0/data/training/annotation_vectors.npy", annotations_vector_list)

#######################################################################################################################

# Extract data from AtTRACT database for Finetuning with RNA-binding Proteins
print("Extracting RBP-Sequences and RNA-Motifs starting.")

# Load AtTRACT database
attract_data = pd.read_csv(r"./data/attract/ATtRACT_db.txt", sep = "\t")

# Iterating through every line of AtTRACT database to extract binding rna_motif, score and find the corresponding amino acid sequence
annotations_vector_list = []
rbp_rna_motifs_list = []
mask_rbp_rna_motifs_list = []
rbp_sequences_list = []
mask_rbp_sequences_list = []
protein_family_list = []
score_list = []
unified_record_sequence_length = 512
unified_rna_motif_length = 16
id_entry_before = 0

for line in attract_data.iloc():

    # Extraction of rna_motif, protein family and score from ATtract database. Function within the function below finds
    # the corresponding amino acid sequence.
    acc_uniprot, rna_motif, sequence_list, mask_rbp_sequence_list, taxonomy, ncbi_taxid, go_terms, family, score, id_entry_before = att.extract_data_from_attract_database(line, record_dict, unified_record_sequence_length, unified_rna_motif_length, id_entry_before)

    rna_motif = str(rna_motif)
    rna_motif = sequ.translate_rna_motif_in_anno_vec(rna_motif)

    if acc_uniprot != 0:

        annotations_vector = np.array(annotations_vector_saved, dtype=float)

        counter = 0

        for go_term in go_terms:
            if go_term not in go_terms_catalytic_activity_list:
                counter += 1

        if counter == len(go_terms):

            for sequence, mask_rbp_sequence in zip(sequence_list, mask_rbp_sequence_list):

                aa_distribution, distribution_code = anno.get_aa_distribution(sequence)
                aa_types, aa_types_description, total_charge = anno.get_aa_types(sequence)

                # Modification of blank annotations vector based on annotations of the record extracted from the record_dict
                annotations_vector = anno_vec.entry_specific_modification_annotation_vector(acc_uniprot, annotations_vector, annotations_vector_description, go_terms, go_terms_semantic_similarity_dict, taxonomy, total_charge, aa_distribution, distribution_code, aa_types, aa_types_description, bindingdb_uniprot_ids, rna_motif)

                # Sequence is converted into a multi-hot-encoding matrice to be readable for the neural network.
                rows_amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '0']
                #sequence = anno.generation_sequence_matrice(sequence, rows_amino_acids)
                sequence = sequ.generation_sequence_one_hot_preparation(sequence, rows_amino_acids)
                #rows_nucleic_bases = ['A', 'C', 'G', 'U', '0']
                #rna_motif_matrice = anno.generation_sequence_matrice(rna_motif, rows_nucleic_bases)
                #rna_motif_matrice = sequ.generation_sequence_one_hot_preparation(rna_motif[0], rows_nucleic_bases)

                if len(sequence) == unified_record_sequence_length: # and len(rna_motif_matrice) == unified_rna_motif_length:
                    annotations_vector_list.append(annotations_vector)
                    score_list.append(score)
                    #rbp_rna_motifs_list.append(rna_motif_matrice)
                    #mask_rbp_rna_motifs_list.append(mask_rbp_rna_motif)
                    rbp_sequences_list.append(sequence)
                    mask_rbp_sequences_list.append(mask_rbp_sequence)
                    protein_family_list.append(family)

                annotations_vector = np.array(annotations_vector_saved, dtype=float)

# Sequences, rna motifs and scores are sorted for the protein families.
protein_family_list_argsort = np.argsort(protein_family_list)
annotations_vector_list = anno.sorting_entries(protein_family_list_argsort, annotations_vector_list)
protein_family_list = anno.sorting_entries(protein_family_list_argsort, protein_family_list)
rbp_sequences_list = anno.sorting_entries(protein_family_list_argsort, rbp_sequences_list)
mask_rbp_sequences_list = anno.sorting_entries(protein_family_list_argsort, mask_rbp_sequences_list)
#rbp_rna_motifs_list = anno.sorting_entries(protein_family_list_argsort, rbp_rna_motifs_list)
#mask_rbp_rna_motifs_list = anno.sorting_entries(protein_family_list_argsort, mask_rbp_rna_motifs_list)
score_list = anno.sorting_entries(protein_family_list_argsort, score_list)

annotations_vector_list = np.asarray(annotations_vector_list)
protein_family_list = np.asarray(protein_family_list)
#rbp_rna_motifs_list = np.asarray(rbp_rna_motifs_list)
#mask_rbp_rna_motifs_list = np.asarray(mask_rbp_rna_motifs_list)
rbp_sequences_list = np.asarray(rbp_sequences_list)
mask_rbp_sequences_list = np.asarray(mask_rbp_sequences_list)
score_list = np.asarray(score_list)

# Saving the rna_motif, binding sequences and scores for use in the neural network
np.save("/beegfs/work/ws/hd_vu199-preprocessing_arnoldt-0/data/finetuning/rbp_annotation_vectors.npy", annotations_vector_list)
np.save("/beegfs/work/ws/hd_vu199-preprocessing_arnoldt-0/data/finetuning/protein_families.npy", protein_family_list)
#np.save("/beegfs/work/ws/hd_vu199-preprocessing_arnoldt-0/data/finetuning/rbp_rna_motifs.npy", rbp_rna_motifs_list)
#np.save("/beegfs/work/ws/hd_vu199-preprocessing_arnoldt-0/data/finetuning/mask_rbp_rna_motifs.npy", mask_rbp_rna_motifs_list)
np.save("/beegfs/work/ws/hd_vu199-preprocessing_arnoldt-0/data/finetuning/rbp_sequences.npy", rbp_sequences_list)
np.save("/beegfs/work/ws/hd_vu199-preprocessing_arnoldt-0/data/finetuning/mask_rbp_sequences.npy", mask_rbp_sequences_list)
np.save("/beegfs/work/ws/hd_vu199-preprocessing_arnoldt-0/data/finetuning/scores.npy", score_list)
