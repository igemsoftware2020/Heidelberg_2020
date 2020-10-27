# PRISM

## Model

The state dict of pretrained neural network model "state_dict_pretraining.pt" as well as the state dict of the finetuned neural network model "state_dict_finetuning.pt" can be accessed via the folder "data". For more information concerning the model architecture please consult the [iGEM Heidelberg 2020 Wiki page](https://2020.igem.org/Team:Heidelberg/Software/PRISM).

## Initial training
Calling the file ```pretraining.py``` trains the Neural Network on all available entries from the Swissprot database and further mentioned databases. For this the encoder is presented the sequences, masks and annotation_vectors generated in the preprocessing pipeline. The trained model’s state dict is automatically saved.

## Finetuning 
Calling the file ```finetuning.py``` lets the Neural Network enters the finetuning step, where the pretrained model’s state dict is automatically loaded. For Finetuning all entries with their rna-motifs in Attract Database are used if the entry is present in Swissprot database. The trained model’s state dict is automatically saved.
