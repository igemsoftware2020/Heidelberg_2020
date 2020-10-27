import sys
import random

from matplotlib import pyplot as plt

import torch
import torch.nn as nn
import torch.nn.functional as func
from torchsupport.training.training import SupervisedTraining
from torchsupport.structured import DataParallel as SDP
from torchsupport.structured import PackedTensor

from rnasupport.data.refold import RNAGraphFoldData, RNAGraphFoldRecomputeData
from rnasupport.modules.baseline import ConvAll, ConvBonds

def valid_callback(training, data, predictions):
  inputs, (labels, mask) = data
  confusion = torch.zeros(4, 4)
  labels = labels[mask]
  predictions = predictions[0][0][mask]
  for label, prediction in zip(labels, predictions):
    pred = prediction.argmax(dim=0)
    confusion[label, pred] += 1
  fig, ax = plt.subplots()
  ax.imshow(confusion / confusion.sum(dim=1, keepdim=True), cmap="Reds")
  training.writer.add_figure("confusion", fig, training.step_id)

class DebugLoss(nn.CrossEntropyLoss):
  def forward(self, inputs, target):
    target, mask = target
    return super().forward(inputs[mask], target[mask])

class DDP(nn.Module):
  def __init__(self, net):
    super().__init__()
    self.net = net

  def forward(self, *args):
    inputs = []
    for arg in args:
      if isinstance(arg, PackedTensor):
        inputs.append(arg.tensor)
      else:
        inputs.append(arg)
    return self.net(*inputs)

if __name__ == "__main__":
  #with torch.autograd.detect_anomaly():
    data = RNAGraphFoldRecomputeData()
    valid_data = RNAGraphFoldRecomputeData()
    net = SDP(ConvBonds())
    training = SupervisedTraining(
      net, data, valid_data,
      [DebugLoss()],
      batch_size=32,
      max_epochs=1000,
      optimizer=lambda x: torch.optim.Adam(x, lr=1e-4),#, weight_decay=1e-4),
      device="cuda:0",
      network_name="baseline/infinite-data-poe-2",
      valid_callback=valid_callback
    )
    final_net = training.train()
