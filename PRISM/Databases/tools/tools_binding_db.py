def extraction_uniprot_ids_bindingdb():
  """
  Function to open BindingDBUniprot.txt and extracting the uniprot ids for proteins which bind ligands.

  :return bindingdb_uniprot_ids: uniprot accession ids for proteins which bind ligands
  """
  bindingdb_uniprot_ids = []

  with open("./Databases/data/binding_db/BindingDB_UniProt.txt", "r") as bindingdb_uniprot:
    for line in bindingdb_uniprot:
      stripped_line = line.strip()
      print(stripped_line)

      stripped_line = stripped_line.replace(",", " ")
      stripped_line = stripped_line.replace(", ", " ")
      stripped_line = stripped_line.replace(" ,", " ")
      stripped_line = stripped_line.replace(" , ", " ")

      _, _, uniprotids = stripped_line.partition("\t")

      for element in range(0, len(uniprotids.partition(" "))):
        if uniprotids.partition(" ")[element] != "" and uniprotids.partition(" ")[element] != " ":
          bindingdb_uniprot_ids.append(uniprotids.partition(" ")[element])

  bindingdb_uniprot.close()

  bindingdb_uniprot_ids = bindingdb_uniprot_ids[1:]

  return bindingdb_uniprot_ids

def update_ligand_binding(id, bindingdb_uniprot_ids):
  """
  Function checks with a protein's id whether the protein is binding a ligand by comparing with bindingdb_uniprot_ids.

  :param id: uniprot accession number of the handled record
  :param bindingdb_uniprot_ids: list of uniprot accession numbers which bind a ligand
  :return binding_db_bool: boolean, whether id is in bindingdb_uniprot_ids
  """

  binding_db_bool = False

  if id in bindingdb_uniprot_ids:
    binding_db_bool = True

  return binding_db_bool
