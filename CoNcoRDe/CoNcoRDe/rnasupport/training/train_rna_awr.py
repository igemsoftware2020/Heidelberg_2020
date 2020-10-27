import sys
import random

from matplotlib import pyplot as plt

import torch
import torch.nn as nn
import torch.nn.functional as func
from torchsupport.training.training import SupervisedTraining
from torchsupport.structured import DataParallel as SDP
from torchsupport.interacting.shared_data import SharedModule
from torchsupport.interacting.policies.basic import CategoricalPolicy
from torchsupport.interacting.crr import CRRTraining
from torchsupport.interacting.bdpi import BDPITraining

from rnasupport.environments.block_rna import SeqRNA
from rnasupport.modules.policy import ConnectedConvPolicy

if __name__ == "__main__":
  policy = ConnectedConvPolicy(64, depth=2)
  value = ConnectedConvPolicy(64, depth=2)
  agent = CategoricalPolicy(SharedModule(policy))
  env = SeqRNA(min_size=32, max_size=33, close=True)
  training = CRRTraining(
    policy, value, agent, env,
    network_name=f"awrna/test-fixed-5",
    verbose=True, beta=1.0,
    auxiliary_steps=1,
    discount=0.990,
    clip=20,
    n_workers=8,
    device="cpu",
    batch_size=32,
    double=True,
    buffer_size=20_000,
    policy_steps=1
  )

  # training = BDPITraining(
  #   policy, value, agent, env,
  #   network_name=f"awrna/bdpi-1",
  #   discount=0.99,
  #   clones=4,
  #   critic_updates=4,
  #   gradient_updates=20,
  #   batch_size=32,
  #   double=True,
  #   buffer_size=20_000,
  #   device="cuda:0",
  #   verbose=True
  # )

  training.train()
