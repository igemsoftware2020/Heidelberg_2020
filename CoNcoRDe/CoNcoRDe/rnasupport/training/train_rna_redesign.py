import sys
import random

from matplotlib import pyplot as plt

import torch
import torch.nn as nn
import torch.nn.functional as func
from torchsupport.training.training import SupervisedTraining
from torchsupport.structured import DataParallel as SDP

from rnasupport.data.refold import RNAGraphRefoldData
from rnasupport.modules.baseline import RedesignConvBond

def valid_callback(training, data, predictions):
  inputs, labels = data
  confusion = torch.zeros(4, 4)
  for label, prediction in zip(labels, predictions[0][0]):
    pred = prediction.argmax(dim=0)
    confusion[label, pred] += 1
  fig, ax = plt.subplots()
  ax.imshow(confusion / confusion.sum(dim=1, keepdim=True), cmap="Reds")
  training.writer.add_figure("confusion", fig, training.step_id)

class DebugLoss(nn.CrossEntropyLoss):
  def forward(self, inputs, target):
    print(inputs, target)
    return super().forward(inputs, target)

if __name__ == "__main__":
  #with torch.autograd.detect_anomaly():
    data = RNAGraphRefoldData(sys.argv[1], mode="Train")
    valid_data = RNAGraphRefoldData(sys.argv[1], mode="Valid")
    net = SDP(RedesignConvBond(depth=1))
    training = SupervisedTraining(
      net, data, valid_data,
      [DebugLoss()],
      batch_size=32,
      max_epochs=1000,
      optimizer=lambda x: torch.optim.Adam(x, lr=1e-3, weight_decay=1e-4),
      device="cuda:0",
      network_name="baseline/refold-4",
      valid_callback=valid_callback
    ).load()
    final_net = training.train()
