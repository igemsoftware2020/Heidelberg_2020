import torch
import torch.nn as nn

def loss_mask(criterion, pred, target, mask):
  
  idx = (mask==0).view(-1).nonzero().view(-1)
  
  pred = pred.transpose(1,2).reshape(-1,21)
  pred = pred[idx]
  target = target.view(-1)[idx]
  loss = criterion(pred, target)
  
  return loss
