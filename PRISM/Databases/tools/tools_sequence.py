from Bio import SeqIO
import random
import numpy as np

def extract_sequence_from_fasta(acc_uniprot):
   """
   Function extracts the sequence from a fasta file for a specific uniprot accession number.

   :param acc_uniprot: the uniprot accession number of the record which is dealt with in the following
   :return sequence: the sequence of the handled record
   """

   for record in SeqIO.parse("./data/attract/" + acc_uniprot + ".fasta", "fasta"):
      sequence = str(record.seq)

   return sequence


def unify_sequence_length(sequence, unified_record_sequence_length):
   """
   Function takes the sequence of the handled record and adds amino acids until the unified length of all sequences in the neural network is reached.

   :param sequence: the sequence of the handled record
   :param unified_sequence_length: the unified length of all sequences used in the neural network
   :return sequence: the sequence of the handled record with a length of the unified_record_sequence_length
   """

   missing_sequence_length = unified_record_sequence_length - len(sequence)

   sequence = sequence + missing_sequence_length * "0"

   return sequence


def cut_sequence_length(sequence, max_sequence_length):
   """
   Function takes the sequence of the handled record and cuts the sequence when the max length of all sequences used in the neural network is reached.

   :param sequence: the sequence of the handled record
   :param max_sequence_length: max length of all sequences used in the neural network
   :return sequence_list: a list with the sequences with a length of the max sequence_length
   :return mask_list: a list of masks for the sequences in the sequence_list
   """

   multiplier = divmod(len(sequence), max_sequence_length)

   sequence_list = []
   mask_list = []

   if multiplier[0] > 0:
      multiplier_loop = multiplier[0]
      rest = 0
      if multiplier[1] > 0:
         multiplier_loop = multiplier_loop + 1
         for i in range(0, multiplier_loop):
            if i == multiplier_loop:
               rest = multiplier[1]
            sequence_list.append(sequence[0+max_sequence_length*i:max_sequence_length+max_sequence_length*i+rest])
            mask_list.append(sequence_mask(sequence[0+max_sequence_length*i:max_sequence_length+max_sequence_length*i+rest], max_sequence_length))
      else:
         sequence_list.append(sequence)
         mask_list.append(sequence_mask(sequence, max_sequence_length))
   else:
      sequence_list.append(sequence)
      mask_list.append(sequence_mask(sequence, max_sequence_length))

   if len(sequence_list[-1]) < max_sequence_length:
      sequence_list[-1] = unify_sequence_length(sequence_list[-1], max_sequence_length)

   return sequence_list, mask_list

def generation_sequence_one_hot_preparation(sequence, rows):
   '''
   Function converts the amino acids in a sequence in a number corresponding to their position in rows.

   :param sequence: the sequence of the handled record
   :param rows: string of all amino acids used in the pipeline.
   :return sequence_preparation: sequence with numbers representing an amino acid.
   '''

   sequence_preparation = []

   for element in sequence:
      sequence_preparation.append(rows.index(element))

   return sequence_preparation


def sequence_mask(sequence, unified_record_sequence_length):
   '''
   Function creates sequence mask for Neural Network to mask padded sequence elements.

   :param sequence: the sequence of the handled record
   :param unified_record_sequence_length: the unified length of all sequences used in the neural network
   :return sequence_mask: sequence mask for the sequence of the handled record
   '''

   sequence_mask = [0] * len(sequence)
   sequence_mask_padding = [1] * (unified_record_sequence_length - len(sequence))

   sequence_mask = sequence_mask + sequence_mask_padding

   return sequence_mask


def resolve_amino_acids(sequence):
   '''
   Function replaces the amino acids 'O', 'U', 'B', 'Z', 'X' and 'J' either with their natural amino acid counterpart or a random choice of all possible amino acids.

   :param sequence: the sequence of the handled record
   :return sequence: the modified sequence of the handled record
   '''

   sequence = list(sequence)

   for element in range(0, len(sequence)):
      if sequence[element] == 'O':
         sequence[element] = 'K'
      elif sequence[element] == 'U':
         sequence[element] = 'C'
      elif sequence[element] == 'B':
         sequence[element] = random.choice(['N', 'D'])
      elif sequence[element] == 'Z':
         sequence[element] = random.choice(['Q', 'E'])
      elif sequence[element] == 'X':
         sequence[element] = random.choice(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'])
      elif sequence[element] == 'J':
         sequence[element] = random.choice(['L', 'I'])

   sequence = "".join(sequence)

   return sequence


def translate_rna_motif_in_anno_vec(rna_motif):
   '''
   Function translates an RNA-Motif to rna_motif_anno_vec as a one-hot encoded lsit, which can be included in the annotations_vector.

   :param rna_motif: RNA-motif bound by a RNA-binding protein.
   :return rna_motif_anno_vec: list of RNA-motif prepared to be included in a annotations_vector.
   '''

   rna_motif_anno_vec = [0.0] * 4 * 12
   
   base_shift = 0
   
   if rna_motif != "":
      for element in range(0, len(rna_motif)):

         base = rna_motif[element]

         if base == "A":
            base_shift = 0
         elif base == "C":
            base_shift = 1
         elif base == "G":
            base_shift = 2
         elif base == "U":
            base_shift = 3

         rna_motif_anno_vec[0 + 4 * element + base_shift] = 1.0

         base_shift = 0
         
   rna_motif_anno_vec = np.asarray(rna_motif_anno_vec, dtype=float)

   return rna_motif_anno_vec
