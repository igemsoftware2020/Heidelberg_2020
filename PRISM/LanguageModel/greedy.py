import torch
import torch.nn.functional
import numpy as np
import time



# Dictionaries
as_to_idx = {0:"A", 
             1:"R", 
             2:"N", 
             3:"D", 
             4:"C", 
             5:"Q", 
             6:"E", 
             7:"G", 
             8:"H", 
             9:"I", 
             10:"L", 
             11:"K", 
             12:"M", 
             13:"F", 
             14:"P", 
             15:"S", 
             16:"T", 
             17:"W", 
             18:"Y",
             19:"V",
             20:"padding"}
as_to_idx2 = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "padding"]
as_to_idx3 = 'ARNDCQEGHILKMFPSTWYV0'


# For one sequence [1, 21, 512]
def greedy(logits, temperature, dictionary):
    """
    It samples amino acids from a prediction output, given your vocabulary/ idx to word definition and a temperature.
    It first computes the temperature based softmax and then uses argmax to compute the index with the highest probability.

    Args:
        logits: prediction output of your model.
        temperature: for calculating the temperature based softmax.
        dictionary: index to word definition
    """

    logits = torch.squeeze(logits)
    logits = logits.transpose(0, 1)
    # Remove first aminoacids = padding
    # logits = logits[:, :19]

    logits = torch.nn.functional.softmax(logits/temperature, dim = 1)

    seq = torch.argmax(logits, dim=1)
    
    seq = seq.cpu().numpy()

    translation = []
    for amino in  seq:
        
        amino = dictionary[amino]
        translation.append(amino)

        print(amino, sep="", end="", flush=True)
        time.sleep(0.01)
    return translation

        
    
# Create a bag of generated proteins [10, 21, 512]
def greedy_bag(logits, temperature, dictionary):
    """
    Given a tensor of x predictions for x sequences. The function returns a bag of generated proteins.

    Args:
        logits: prediction output of your model.
        temperature: for calculating the temperature based softmax.
        dictionary: index to word definition
    """
    logits = logits.transpose(1, 2)

    # Remove padding
    # logits = logits[:, :, :19]

    bag_of_proteins ={}
    idx = 0

    for prob in logits:

         prob = torch.nn.functional.softmax(prob/temperature, dim = 1)

         seq = torch.argmax(prob, dim=1)
    
         seq = seq.cpu().numpy()

         translation = []

         for amino in  seq:
        
             amino = dictionary[amino]
             translation.append(amino)
             # print(amino, sep="", end="", flush=True)
             # time.sleep(0.01)
            
         bag_of_proteins[f'{idx}'] = translation
         idx +=1

    return bag_of_proteins

# Translate given indices
def amino_code_bag(indices, dictionary):
    """
    Given a sequence of indices in the range of your dictionary, the functions translates the indices to the matching tokken.
    
    Args:
        indices = list of indices
        dictionary = index to word definition
    """
    indices = indices.cpu().numpy()
    bag_of_proteins = {}
    idx = 0

    for protein in indices:
    
         translation = [] 

         for amino in protein:
    
            amino = dictionary[amino]
            translation.append(amino)

         bag_of_proteins[f'{idx}'] = translation 
         idx +=1
         
    return bag_of_proteins


"""
def nucleus_sampling(logits, top_p):
        # source: https://arxiv.org/pdf/1904.09751.pdf
        # source: https://gist.github.com/thomwolf/1a5a29f6962089e871b94cbd09daf317
    if top_p > 0.0:
        sorted_logits, sorted_indices = torch.sort(logits, descending=True)
        cumulative_probs = torch.cumsum(F.softmax(sorted_logits, dim=-1), dim=-1)

        # Remove tokens with cumulative probability above the threshold
        sorted_indices_to_remove = cumulative_probs > top_p
        # Shift the indices to the right to keep also the first token above the threshold
        sorted_indices_to_remove[..., 1:] = sorted_indices_to_remove[..., :-1].clone()
        sorted_indices_to_remove[..., 0] = 0

        indices_to_remove = sorted_indices[sorted_indices_to_remove]
        logits[indices_to_remove] = filter_value
    return logits
"""