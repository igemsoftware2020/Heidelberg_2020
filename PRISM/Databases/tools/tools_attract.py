from tools import tools_uniprot as uni
from tools import tools_sequence as sequ
from tools import tools_swissprot as swiss
import pandas as pd

def extract_data_from_attract_database(line, record_dict, unified_record_sequence_length, unified_rna_motif_length, id_entry_before):
   """
   Function extracts annotations from attract database file.

   :param line: line of csv file of attract database file containing information for one record
   :param record_dict: dictionary of all records of swissprot
   :param unified_record_sequence_length: the unified length of all sequences used in the neural network
   :param unified_rna_motif_length: the unified length of all rna motifs used in the neural network
   :param id_entry_before: id of the entry upfront of this entry
   :return acc_uniprot: uniprot accession number of record
   :return rna_motif: binding motif of record
   :return sequence_list: sequences of record
   :return mask_rbp_sequences_list: mask for sequences of record for the neural network
   :return taxonomy: taxonomy of record
   :return ncbi_taxid: ncbi_taxid of record
   :return go_terms: go terms of record
   :return family: protein family of record
   :return score: score describing binding affinity of RNA-motif
   """
   
   id = line["Gene_id"]
   id = str(id)
   id_entry_attract = id
   organism = line["Organism"]
   rna_motif = line["Motif"]
   family = line["Family"]
   data_source = line["Database"]
   score = line["Score"]
   if score[-2:] == "**":
      score = score[:-2]
   score = float(score)

   if id != id_entry_before:
      
      # Calling function, which gives back the right database identifier, which can be used for uniprot id mapping
      database, id = identify_database(id, data_source)

      if database != 0:

         #EXAMPLE
         #database = "ENSEMBL_ID"
         #id = "ENSG00000148584"
         #DELETE

         # Calling function, which gives back the NCBI taxonomy identifier (number) for a specific organism.
         ncbi_taxid = uni.convert_organism_to_ncbi_taxid(organism)

         # Calling function, which accesses uniport id mapping service. This function then outputs the uniprot identifier for a identifier of another database
         record_list = uni.uniprot_id_mapping(database, id)

         if record_list != 0:

            # Calling function, which determines the right uniprot entry based on the NCBI taxonomy identifier
            acc_uniprot, sequence_list, mask_rbp_sequences_list, taxonomy, ncbi_taxid, go_terms, family = determine_right_uniprot_entry(record_list, ncbi_taxid, record_dict, unified_record_sequence_length)

            if acc_uniprot != 0:
               # rna_motif, mask_rna_motif = sequ.cut_sequence_length(rna_motif, unified_rna_motif_length)

               return acc_uniprot, rna_motif, sequence_list, mask_rbp_sequences_list, taxonomy, ncbi_taxid, go_terms, family, score, id_entry_attract

            else:
               return 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

         else:
            return 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

      else:
         return 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
   
   else:
      return 0, 0, 0, 0, 0, 0, 0, 0, 0, 0


def identify_database(id, data_source):
   """
   Function outputs the database, in which the accession code provided with attract_data can be found in.

   :param id: id of entry in ATtract database
   :param data_source: source of data mentioned in ATtract database
   :return database: database, in which the accession codes ca be found
   :return id: id of entry in ATtract database, in some cases modified.
   """
   #Source: https://ebi12.uniprot.org/help/api_idmapping

   database = 0

   if id[0:3] == "ENS":
      database = "ENSEMBL_ID"
   elif id[0:2] == "FB":
      database = "FLYBASE_ID"
   elif id[0:2] == "WB":
      database = "WORMBASE_ID"
   elif id[0:2] == "AT":
      database = "ARAPORT_ID"
   elif id[0:1] == "B":
      database = "STRING_ID"
      id = "7091." + id + "-TA"
   elif id[0:2] == "XB":
      database = "XENBASE_ID"
   elif id[0:2] == "YO":
      database = "GENENAME"
   elif id[0:10] == "NEMVEDARFT":
      database = "GENENAME"
      id = id[11:]
   elif id[0:2] == "GR":
      database = "MAIZEGDB_ID"
   elif id[0:3] == "SMP":
      database = "GENENAME"
   elif id[0:3] == "TVAG":
      database = "GENENAME"
   elif data_source == "PDB": #For "M", "D", and "GL"
      database = "EMBL_ID"

   return database, id


def determine_right_uniprot_entry(record_list, ncbi_taxid, record_dict, unified_record_sequence_length):
   """
   Function extracts record of list of records (result of uniprot id mapping service) based on the ncbi taxonomy id.

   :param record_list: list of records (result of uniprot id mapping service)
   :param ncbi_taxid: ncbi taxid for record
   :param record_dict: dictionary of all records of swissprot
   :param unified_record_sequence_length:
   :return acc_uniprot: uniprot accession number of record
   :return sequence: sequence of record
   :return taxonomy: taxonomy of record
   :return ncbi_taxid: ncbi_taxid of record
   :return go_terms: go terms of record
   """

   sequence_list = 0
   mask_list_sequence_list = 0
   taxonomy = 0
   go_terms = 0
   family = 0

   # TESTING - DELETE LATER
   counter1 = 0
   counter2 = 0

   for acc_uniprot in record_list:
      try:
         record = record_dict[acc_uniprot]

         # TESTING - DELETE LATER
         counter1 = counter1 + 1

         if record.annotations["ncbi_taxid"][0] == ncbi_taxid:

            _, sequence_list, mask_list_sequence_list, taxonomy, _, go_terms, family = swiss.extraction_annotations(record_dict, acc_uniprot, unified_record_sequence_length)

            # TESTING - DELETE LATER
            counter2 = counter2 + 1

         # TESTING - DELETE LATER
         if counter1 > 1:
            print("To many reviewed ones. ANOTHER CRITERIUM NECESSARY!")
         if counter2 > 1:
            print("To many right ones. ANOTHER CRITERIUM NECESSARY!")

      except:
         continue

   if sequence_list == 0:
      acc_uniprot = 0

   return acc_uniprot, sequence_list, mask_list_sequence_list, taxonomy, ncbi_taxid, go_terms, family
