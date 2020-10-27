import os
import sys
sys.path.append(os.getcwd())
sys.path.append("./LanguageModel")

import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np

import sets
import preprocessing_pca as pca
import dataset
import greedy
import rezero_transformer

def logo():
    print('*********************************************************************')
    print('\
*                _____  _____  _____  _____ __  __                  *\n\
*               |  __ \|  __ \|_   _|/ ____|  \/  |                 *\n\
*               | |__) | |__) | | | | (___ | \  / |                 *\n\
*               |  ___/|  _  /  | |  \___ \| |\/| |                 *\n\
*               | |    | | \ \ _| |_ ____) | |  | |                 *\n\
*               |_|    |_|  \_\_____|_____/|_|  |_|                 *')
    print('*                                                                   *')
    print("* Authors: Wangjun, Hu; Arnoldt, Lucas; Chernova, Elizaveta         *")
    print("* iGEM Team Heidelberg 2020                                         *")
    print("* Please email your comments to: igemhd@protonmail.com              *")
    print('*********************************************************************')

def generation(annotations_vector):
    
    logo()

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    encoder = rezero_transformer.Encoder(d_model_size=512)
    encoder = nn.DataParallel(encoder)
    PATH = "./LanguageModel/data/state_dict_finetuning.pt"
    encoder.load_state_dict(torch.load(PATH))
    encoder.to(device)

    criterion = torch.nn.CrossEntropyLoss()
    encoder.eval()

    annotations_vector = annotations_vector
    sequence = np.array([[0.0] * 512])
    
    mask = np.array([[0.0] * 512])
    dataset_generation = dataset.UniprotDataset(sequence, annotations_vector, mask)
    
    for idx, data in enumerate(dataset_generation):
        sequence, annotations_vector, mask = data
        
        sequence, annotations_vector, mask = sequence.to(device), annotations_vector.to(device), mask.to(device)
        sequence = torch.unsqueeze(sequence, 0)
        annotations_vector = torch.unsqueeze(annotations_vector, 0)
        mask = torch.unsqueeze(mask, 0)
        print(sequence.shape)
        print(annotations_vector.shape)
        print(mask.shape)
        output = encoder(seq = sequence, anno_vec = annotations_vector, mask = mask)

        # generation
        as_to_idx = 'ARNDCQEGHILKMFPSTWYV0'
        generated_seq = greedy.greedy(logits = output, temperature = 0.7, dictionary = as_to_idx)
    
    return generated_seq
