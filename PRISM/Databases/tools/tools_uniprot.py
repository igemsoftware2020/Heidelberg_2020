import urllib.parse
import urllib.request
import requests

from tools import tools_sequence as sequ

def uniprot_id_mapping(database, accession_code):
   """
   Function provides programmatic access to Uniprot ID Mapping. It converts the provided accession_code for a certain database to a uniprot accession code, which is returned.

   :param database: database where the accession code can be found in.
   :param accession_code: accession code for the database
   :return acc_uniprot: uniprot accession code
   """

   print(database)
   print(accession_code)   

   url = 'https://ebi12.uniprot.org/uploadlists/' #OLD: 'https://www.uniprot.org/uploadlists/'

   params = {
      'from': database,
      'to': 'ACC',
      'format': 'tab',
      'query': accession_code,
      'status': 'reviewed'
   }

   data = urllib.parse.urlencode(params)
   data = data.encode('utf-8')
   req = urllib.request.Request(url, data)
   with urllib.request.urlopen(req) as f:
      response = f.read()

   response = str(response).split("b'From\\tTo")
   response = response[1].split("\\n'")
   response = response[0].split("\\n" + accession_code + "\\t")
   response_list = response [1:]

   if len(response_list) > 0:
      return response_list
   else:
      return 0

def uniprot_download_single_entry(acc_uniprot):
   """
   Function downloads the corresponding .fasta file for a uniprot accession number and saves it in the folder ./data/rbps.

   :param acc_uniprot: uniprot accession number
   :return sequence_list: list of sequences
   """

   url_entry_uniprot = "https://www.uniprot.org/uniprot/" + acc_uniprot + ".fasta"
   request_entry_uniprot = requests.get(url_entry_uniprot, allow_redirects=True)
   open("./data/attract/" + acc_uniprot + ".fasta", 'wb').write(request_entry_uniprot.content)

   sequence = sequ.extract_sequence_from_fasta(acc_uniprot)

   max_sequence_length = 500
   sequence_list = sequ.cut_sequence_length(sequence, max_sequence_length)

   return sequence_list

def convert_organism_to_ncbi_taxid(organism):
   """
   Function converts an organism input from the AtTRACT database to the NCBI-TaxId.

   :param organism: designation of the organism a protein (record) is found in
   :return ncbi_taxid: identifier for a specific organism by ncbi
   """

   ncbi_taxid_list = {
      "Mus_musculus": '10090',
      "Homo_sapiens": '9606',
      "Gallus_gallus": '9031',
      "Drosophila_melanogaster": '7227',
      "Rattus_norvegicus": '10116',
      "Spodoptera_frugiperda": '7108',
      "Aspergillus_nidulans": '162425',
      "Caenorhabditis_elegans": '6239',
      "Arabidopsis_thaliana": '3702',
      "Bombyx_mori": '7091',
      "Brachypodium_distachyon": '15368',
      "Xenopus_laevis": '8355',
      "Chaetomium_thermophilum": '209285',
      "Phytophthora_ramorum": '164328',
      "Saccharomyces_cerevisiae": '4932',
      "Leishmania_major": '5664',
      "Naegleria_gruberi": '5762',
      "Physcomitrella_patens": '3218',
      "Ostreococcus_tauri": '70448',
      "Oryctolagus_cuniculus": '9986',
      "Tetraodon_nigroviridis": '99883',
      "Xenopus_tropicalis": '8364',
      "Vanderwaltozyma_polyspora": '36033',
      "Plasmodium_falciparum": '5833',
      "Saccharomyces_cerevisiae_s288c": '1196866',
      "Mesocricetus_auratus": '10036',
      "Cricetulus_griseus": '10029',
      "Neurospora_crassa": '5141',
      "Nematostella_vectensis": '45351',
      "Danio_rerio": '7955',
      "Zea_mays": '4577',
      "Bos_taurus": '9913',
      "Rhizopus_oryzae": '64495',
      "Schistosoma_mansoni": '6183',
      "Trypanosoma_brucei": '5691',
      "Oryzias_latipes": '8090',
      "Thalassiosira_pseudonana": '35128',
      "Trichomonas_vaginalis": '5722'
   }

   ncbi_taxid = ncbi_taxid_list[organism]

   return ncbi_taxid

