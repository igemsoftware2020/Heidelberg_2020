import torch

from torchsupport.structured import ScatterStructure

class DistanceStructure(ScatterStructure):
  def message(self, source, target):
    difference = (self.connections.float() - self.indices.float()) / 1000
    source_sin = difference.sin()[:, None]
    source_cos = difference.cos()[:, None]
    target_sin = (-difference).sin()[:, None]
    target_cos = (-difference).cos()[:, None]
    source = torch.cat((source[self.connections], source_sin, source_cos), dim=1)
    #target = torch.cat((target[self.indices], target_sin, target_cos), dim=1)
    target = target[self.indices]
    return source, target, self.indices, self.node_count
