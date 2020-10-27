import torch
import torch.nn as nn
import torch.nn.functional as func

from torchsupport.data.namedtuple import namedtuple

from rnasupport.modules.autoregressive import ConnectedConvBlock

class ConnectedConvPolicy(nn.Module):
  data_type = namedtuple("Data", ["logits", "outputs"])
  def __init__(self, size, out_size=4, kernel_size=5, position_size=10, depth=4):
    super().__init__()
    self.position_size = position_size
    self.encode = nn.Conv1d(5, size, 1)
    self.current_blocks = nn.ModuleList([
      ConnectedConvBlock(
        size,
        kernel_size=kernel_size,
        position_size=position_size,
        causal=False
      )
      for idx in range(depth)
    ])
    self.target_blocks = nn.ModuleList([
      ConnectedConvBlock(
        size,
        kernel_size=kernel_size,
        position_size=position_size,
        causal=False
      )
      for idx in range(depth)
    ])
    self.predict = nn.Linear(size, out_size)

  def schema(self):
    return self.data_type(
      logits=torch.zeros(4, dtype=torch.float),
      outputs=None
    )

  def forward(self, state, hidden=None):
    sequence, current, target, mask, place = state
    out = self.encode(sequence)
    for cblock, tblock in zip(self.current_blocks, self.target_blocks):
      out = cblock(out, current, mask)
      out = tblock(out, target, mask)
    ind = torch.arange(sequence.size(0), dtype=torch.long, device=sequence.device)
    out = out.transpose(1, 2)[ind, place.view(-1)]
    out = self.predict(out)
    return out

class PositionPolicy(ConnectedConvPolicy):
  def __init__(self, size, out_size=4, kernel_size=5, position_size=10, depth=4,
               max_size=33):
    super().__init__(
      size, out_size=out_size, kernel_size=kernel_size,
      position_size=position_size, depth=depth
    )
    self.max_size = max_size
    self.predict = nn.Conv1d(size, out_size, 1)

  def schema(self):
    return self.data_type(
      logits=torch.zeros(4 * self.max_size, dtype=torch.float),
      outputs=None
    )

  def forward(self, state, hidden=None):
    sequence, current, target, mask, _ = state
    out = self.encode(sequence)
    for cblock, tblock in zip(self.current_blocks, self.target_blocks):
      out = cblock(out, current, mask)
      out = tblock(out, target, mask)
    out = self.predict(out)
    out[~mask[:, None].expand_as(out)] = -float("inf")
    out = out.view(out.size(0), -1)
    return out
