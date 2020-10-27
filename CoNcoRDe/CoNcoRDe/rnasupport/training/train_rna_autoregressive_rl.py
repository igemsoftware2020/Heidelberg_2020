import sys
import random

from RNA import fold

from matplotlib import pyplot as plt

import torch
import torch.nn as nn
import torch.nn.functional as func
from torchsupport.modules.basic import one_hot_encode
from torchsupport.training.training import SupervisedTraining
from torchsupport.structured import DataParallel as SDP

from rnasupport.data.refold import RNABlockFoldRecomputeData
from rnasupport.modules.autoregressive import ConnectedConvRL
from rnasupport.secstruct import block_bonds

class RLData(RNABlockFoldRecomputeData):
  def __getitem__(self, index):
    size = random.randrange(self.min_size, self.max_size)
    seq = random.choices(["G", "A", "U", "C"], k=size)
    sequence = "".join(seq)
    structure, _ = fold(sequence)

    sequence_gt = torch.zeros(self.max_size, dtype=torch.long)
    sequence_gt[:size] = one_hot_encode(sequence, "GAUC", numeric=True)
    sequence_onehot = torch.zeros(4, self.max_size)
    sequence_onehot[:, :size] = one_hot_encode(sequence, "GAUC")
    structure_gt = torch.zeros(self.max_size, dtype=torch.long)
    structure_gt[:size] = one_hot_encode(structure, ".()", numeric=True)

    mask = torch.zeros(self.max_size, dtype=torch.bool)
    mask[:size] = True
    bonds = torch.zeros(self.max_size, dtype=torch.long)
    causal_bonds = torch.zeros(self.max_size, dtype=torch.long)
    b, cb = block_bonds(structure)
    bonds[:size] = b
    causal_bonds[:size] = cb

    inputs = (
      sequence_onehot, structure_gt, bonds, causal_bonds, mask
    )
    outputs = (
      sequence_gt, mask
    )
    rl = (
      mask
    )

    return inputs, outputs, rl

def valid_callback(training, data, predictions):
  inputs, (labels, mask), _ = data
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

class PGrad(nn.Module):
  def forward(self, inputs, mask):
    logits, rewards = inputs
    result = logits * rewards * mask.float()
    result = result.sum(dim=1).mean(dim=0)
    return result

class DDP(nn.Module):
  def __init__(self, net):
    super().__init__()
    self.net = net

  def forward(self, *args):
    return self.net(*args)

if __name__ == "__main__":
  #with torch.autograd.detect_anomaly():
    data = RLData(max_size=100)
    valid_data = RLData(min_size=10, max_size=100)
    net = SDP(ConnectedConvRL(64, kernel_size=7, depth=4))
    training = SupervisedTraining(
      net, data, valid_data,
      [DebugLoss(), PGrad()],
      batch_size=8,
      max_epochs=1000,
      optimizer=lambda x: torch.optim.Adam(x, lr=1e-4),#, weight_decay=1e-4),
      device="cuda:0",
      network_name="autoregressive-rl/infinite-data-4",
      valid_callback=valid_callback
    )
    final_net = training.train()
