import sys
import random

from matplotlib import pyplot as plt

import torch
import torch.nn as nn
import torch.nn.functional as func
from torchsupport.training.training import SupervisedTraining
from torchsupport.structured import DataParallel as SDP

from rnasupport.data.refold import RNABlockFoldRecomputeData
from rnasupport.modules.autoregressive import ConnectedConvModel, ConnectedConvScheduled

def valid_callback(training, data, predictions):
  inputs, (labels, mask) = data
  confusion = torch.zeros(4, 4)
  labels = labels[mask]
  predictions = predictions[0][0].transpose(0, 1)[:, mask].transpose(0, 1)
  for label, prediction in zip(labels, predictions):
    pred = prediction.argmax(dim=0)
    confusion[label, pred] += 1
  fig, ax = plt.subplots()
  ax.imshow(confusion / confusion.sum(dim=1, keepdim=True), cmap="Reds")
  training.writer.add_figure("confusion", fig, training.step_id)

class DebugLoss(nn.CrossEntropyLoss):
  def forward(self, inputs, target):
    target, mask = target
    return super().forward(
      inputs.transpose(0, 1)[:, mask].transpose(0, 1), target[mask])

class DDP(nn.Module):
  def __init__(self, net):
    super().__init__()
    self.net = net

  def forward(self, *args):
    return self.net(*args)

if __name__ == "__main__":
  #with torch.autograd.detect_anomaly():
    data = RNABlockFoldRecomputeData(max_size=100)
    valid_data = RNABlockFoldRecomputeData(min_size=10, max_size=100)
    net = SDP(ConnectedConvModel(128, kernel_size=9, depth=6))
    training = SupervisedTraining(
      net, data, valid_data,
      [DebugLoss()],
      batch_size=32,
      max_epochs=1000,
      optimizer=lambda x: torch.optim.Adam(x, lr=1e-4),#, weight_decay=1e-4),
      device="cuda:0",
      network_name="autoregressive/infinite-data-17",
      valid_callback=valid_callback
    )
    final_net = training.train()
