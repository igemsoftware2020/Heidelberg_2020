from torch.utils.data import Dataset
import torch
import torch.nn as nn
import numpy as np

class UniprotDataset(Dataset):

    def __init__(self, sequences, annotation_vectors, masks):
        """
        """

        self.sequences = torch.LongTensor(sequences)
        self.annotation_vectors = torch.from_numpy(annotation_vectors)
        self.masks = torch.from_numpy(masks)

        # Create a dictionary mapping each label to a index from 0 to len(classes).
        #self.annotation_vector_to_idx = {x: i for i, x in enumerate(set(self.annotation_vectors))}

    def __len__(self):
        # return length of dataset
        return len(self.sequences)

    def __getitem__(self, idx):

        sequence = self.sequences[idx]
        annotation_vectors = self.annotation_vectors[idx]
        mask = self.masks[idx]

        return sequence, annotation_vectors, mask

'''
class FinetuningDataset(Dataset):

    def __init__(self, sequences, annotation_vectors, rna_motif, masks_sequences, masks_rna_motif):
        """
        """

        self.sequences = torch.LongTensor(sequences)
        self.annotation_vectors = torch.from_numpy(annotation_vectors)
        self.rna_motif = torch.from_numpy(rna_motif)
        self.masks_sequences = torch.from_numpy(masks_sequences)
        self.masks_rna_motif = torch.from_numpy(masks_rna_motif)

    def __len__(self):
        # return length of dataset
        return len(self.sequences)

    def __getitem__(self, idx):
        sequence = self.sequences[idx]
        annotation_vectors = self.annotation_vectors[idx]
        rna_motif = self.arna_motif[idx]
        masks_sequences = self.masks_sequences[idx]
        masks_rna_motif = self.masks_rna_motif[idx]

        return sequence, annotation_vectors, rna_motif, masks_sequences, masks_rna_motif
'''

'''
class LanguageData(Dataset):
  def _getitem_(self, idx):
    sequence = self.data[idx]

    # one hot encode sequence
    sequence_one_hot = one_hot_encode(sequence, AA_CODE)

    # shift sequence
    shifted_sequence = sequence_one_hot.roll(1, dims=1)
    shifted_sequence[:, 0] = 0 # first element should be zero

    membership = SubgraphStructure(torch.zeros(primary.size(0), dtype=torch.long))
    inputs = (PackedTensor(shifted_sequence.permute(1, 0)), membership)
    outputs = PackedTensor(sequence_one_hot.permute(1, 0))

    return inputs, outputs
'''