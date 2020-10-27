import os
import sys
sys.path.append(os.getcwd())
sys.path.append("./data")
sys.path.append("./tools")
sys.path.append("./Databases")
sys.path.append("./Databases/data")
sys.path.append("./LanguageModel")

import torch
import torch.nn as nn
import sets
import numpy as np
import train
from torch.utils.tensorboard import SummaryWriter
import transformer


if __name__ == "__main__":

    set = 5
    multiplier = 2
    size_val_set = 0.1
    
    #Load the datasets prepared in the preprocessing pipeline
    dataset_annotation_vectors = np.load("./Databases/data/finetuning/rbp_annotation_vectors.npy")
    dataset_sequences = np.load("./Databases/data/finetuning/rbp_sequences.npy")
    dataset_masks = np.load("./Databases/data/finetuning/mask_rbp_sequences.npy")
    #dataset_scores = np.load("./Databases/data/finetuning/scores.npy")
    protein_families = np.load("./Databases/data/finetuning/protein_families.npy")
    
    #Determine borders of protein families in datasets
    count_protein_family_list = sets.determine_protein_family_lengths(protein_families)

    #Split the datasets in (1/size_val_set) equally big parts under considerance of the protein families of the entries in the datasets
    dataset_split_annotations_vector = sets.cross_validation_dataset_split(dataset_annotation_vectors, count_protein_family_list, (1 / size_val_set))
    dataset_split_sequences = sets.cross_validation_dataset_split(dataset_sequences, count_protein_family_list, (1 / size_val_set))
    dataset_split_masks = sets.cross_validation_dataset_split(dataset_masks, count_protein_family_list, (1 / size_val_set))
    #dataset_split_scores = sets.cross_validation_dataset_split(dataset_scores, count_protein_family_list, (1 / size_val_set))

    #Initialize the SummaryWriter for tensorboard
    writer = SummaryWriter(f'summary_writer')

    #Initialize the encoder
    encoder = transformer.Encoder(d_model_size=512)
    transformer.freeze_params(encoder)
    transformer.unfreeze_layer(encoder.embedding)
    transformer.unfreeze_layer(encoder.predict)
    encoder = nn.DataParallel(encoder)
    PATH = "./LanguageModel/data/state_dict_pretraining.pt"
    encoder.load_state_dict(torch.load(PATH))
    
    #Set the number of epochs the neural network should be trained
    num_epochs = 20

    train.train(encoder, num_epochs, dataset_split_annotations_vector, dataset_split_sequences, dataset_split_masks, writer, set, multiplier)
