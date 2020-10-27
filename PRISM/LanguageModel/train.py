import os
import sys
sys.path.append(os.getcwd())
sys.path.append("./LanguageModel")

import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import random 

import sets
import preprocessing_pca as pca
import dataset
import loss_mask
import greedy
import confusion_matrix

def train(encoder, num_epochs, dataset_split_annotations_vector, dataset_split_sequences, dataset_split_masks, writer, set, multiplier):
    
    training_set_annotations_vector, validation_set_annotations_vector = sets.cross_validation_dataset_combination(dataset_split_annotations_vector, set)
    training_set_annotations_vector, validation_set_annotations_vector = pca.preprocessing_pca(training_set_annotations_vector, validation_set_annotations_vector)
    training_set_sequences, validation_set_sequences = sets.cross_validation_dataset_combination(dataset_split_sequences, set)
    training_set_masks, validation_set_masks = sets.cross_validation_dataset_combination(dataset_split_masks, set)

    dataset_train = dataset.UniprotDataset(training_set_sequences, training_set_annotations_vector, training_set_masks)
    dataset_val = dataset.UniprotDataset(validation_set_sequences, validation_set_annotations_vector, validation_set_masks)

    batch_size_train = 32
    batch_size_val = batch_size_train

    dataloader_train = torch.utils.data.DataLoader(dataset_train, batch_size=batch_size_train, shuffle=True, num_workers=8, drop_last=True)
    
    #######################################################################################################################

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    encoder.to(device)

    # Loss, optimizer, writer
    criterion = torch.nn.CrossEntropyLoss()
    criterion_val = torch.nn.CrossEntropyLoss(reduction = "sum")
    #criterion = torch.nn.CrossEntropyLoss(weight=torch.from_numpy(training_set_scores))
    step = 0

    # Training
    
    warmup_steps = 50
    d_model = 512
    learning_rate = d_model**(-0.5) * min((1)**(-0.5), 1 * (warmup_steps)**(-1.5))
    optimizer = optim.Adam(encoder.parameters(), lr=learning_rate, betas=(0.9, 0.98))

    optimizer.zero_grad()

    for epoch in range(num_epochs):
         encoder.train()
         print('Epoch {}/{}'.format(epoch, num_epochs - 1))
         print('-' * 10)

         for sequence_batch, annotations_batch, masks_batch in dataloader_train:
             sequence_batch, annotations_batch, masks_batch = sequence_batch.to(device), annotations_batch.to(device), masks_batch.to(device)
             
             target = sequence_batch

             sequence_batch = torch.roll(sequence_batch, 1, dims=0)
             
             sequence_batch[:, 0] = 0
             
             output = encoder(seq = sequence_batch, anno_vec = annotations_batch, mask = masks_batch)

             scores = loss_mask.loss_mask(criterion, output, target, masks_batch)
            
             scores.backward()

             if (step+1) % multiplier == 0:     
                optimizer.step()
                optimizer.zero_grad()
            
             writer.add_scalar('Training loss step', scores, global_step = step)
             
             step +=1
             learning_rate = d_model**(-0.5) * min((step)**(-0.5), step * (warmup_steps)**(-1.5))
             optimizer.param_groups[0]['lr'] = learning_rate
                
             #EVALUATION   
             if step % 100 == 0:
                with torch.no_grad():           
                    encoder.eval()
                    scores = 0
                    dimension = 0
                    y_true = []
                    y_pred = []
                    
                    random_number = random.randint(0,(len(dataset_val)//batch_size_val)-1)
                    random_batch = torch.utils.data.SubsetRandomSampler(list(range(random_number*batch_size_val, (random_number+1)*batch_size_val)))
                    dataloader_val = torch.utils.data.DataLoader(dataset_val, batch_size=batch_size_val, num_workers=8, sampler = random_batch)
                    
                    for sequence_batch, annotations_batch, masks_batch in dataloader_val:
                        print("dataloader_val")
                        sequence_batch, annotations_batch, masks_batch = sequence_batch.to(device), annotations_batch.to(device), masks_batch.to(device)

                        target = sequence_batch
                        sequence_batch = torch.roll(sequence_batch, 1, dims=0)
                        sequence_batch[:, 0] = 0
                        output = encoder(seq = sequence_batch, anno_vec = annotations_batch, mask = masks_batch)
                        
                        tar = target.view(-1)
                        pred = torch.argmax(output, dim=1)
                        pred = pred.view(-1)

                        y_true, y_pred, conf_matrix  = confusion_matrix.create_confusion_matrix(y_true, y_pred, tar, pred)
                        print(conf_matrix)
                        class_names = ["0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"]
                        figure = confusion_matrix.plot_confusion_matrix(conf_matrix, class_names)
                        print(figure)
                        writer.add_figure(f'Confusion matrix-{epoch}', figure, global_step = epoch)

                        scores += loss_mask.loss_mask(criterion_val, output, target, masks_batch)

                        dimension += output.size()[0]

                    writer.add_scalar('Validation loss step', scores/(dimension*output.size()[1]*output.size()[2]), global_step = step)
             
             encoder.train()   

             if step % 1000 == 0:
                torch.save(encoder.state_dict(), f"/beegfs/work/ws/hd_vu199-preprocessing_arnoldt-0/runs/state_dict_encoder_checkpoint-FINETUNING-batchT64_Arnoldt_MASK_WARM50-{epoch}-{step}-{set}.pt")
             
    writer.close()
                 
