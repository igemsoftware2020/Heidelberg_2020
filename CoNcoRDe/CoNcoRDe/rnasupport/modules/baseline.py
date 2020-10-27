import torch
import torch.nn as nn
import torch.nn.functional as func

from torchsupport.structured import NeighbourDotMultiHeadAttention, NeighbourLinear
from torchsupport.structured import scatter

from rnasupport.modules.connections import DistanceStructure

class LocalFeatures(nn.Module):
  def __init__(self, hidden_size=128, depth=4):
    super().__init__()
    self.blocks = nn.ModuleList([
      nn.Conv1d(hidden_size, hidden_size, 3, dilation=idx + 1, padding=idx + 1)
      for idx in range(depth)
    ])
    self.drop = nn.Dropout(0.1)

  def forward(self, inputs):
    out = inputs
    for block in self.blocks:
      out = func.relu(block(out))
    return self.drop(out) + inputs

class NonLocalFeatures(nn.Module):
  def __init__(self, relative_size=2, hidden_size=128, attention_size=128, heads=16):
    super().__init__()
    self.attention = NeighbourLinear(hidden_size + relative_size, hidden_size)
    #self.attention = NeighbourDotMultiHeadAttention(
    #  hidden_size + relative_size, hidden_size, attention_size, heads=heads
    #)
    self.bn = nn.LayerNorm(hidden_size)
    self.drop = nn.Dropout(0.1)

  def forward(self, inputs, structure):
    out = self.bn(inputs)
    out = out
    out = self.attention(out, out, structure)
    return self.drop(out) + inputs

class TransformerBlock(nn.Module):
  def __init__(self, relative_size=2, hidden_size=128, attention_size=128,
               heads=16, depth=4):
    super().__init__()
    self.local = LocalFeatures(hidden_size=hidden_size, depth=depth)
    self.bn = nn.LayerNorm(hidden_size)
    self.connected = NonLocalFeatures(
      relative_size=relative_size, hidden_size=hidden_size,
      attention_size=attention_size, heads=heads
    )

  def forward(self, inputs, rna, structure):
    out = scatter.batched(self.local, self.bn(inputs), rna.indices)
    out = self.connected(out, structure)
    return out

class ConvBonds(nn.Module):
  def __init__(self, in_size=5, depth=4, mlp_depth=4,
               hidden_size=128, attention_size=128,
               relative_size=2, heads=8):
    super().__init__()
    self.blocks = nn.ModuleList([
      TransformerBlock(
        relative_size=relative_size,
        hidden_size=hidden_size,
        attention_size=attention_size,
        heads=heads, depth=mlp_depth
      )
      for idx in range(depth)
    ])

    self.preprocess = nn.Conv1d(in_size + 1, hidden_size, 3, padding=1)
    self.postprocess = nn.Conv1d(hidden_size, 4, 3, padding=1)

  def forward(self, inputs, rna, bonds, stems, unpaired):
    structure = DistanceStructure(0, 0, bonds.indices, bonds.connections, node_count=bonds.node_count)
    pos = torch.arange(inputs.size(0), dtype=torch.float, device=inputs.device) / (1 + torch.repeat_interleave(rna.counts.float(), rna.counts))
    inputs = torch.cat((inputs, pos[:, None]), dim=1)
    out = scatter.batched(self.preprocess, inputs, rna.indices)
    for block in self.blocks:
      out = block(out, rna, structure)
    out = scatter.batched(self.postprocess, out, rna.indices)
    return out

class RedesignConvBond(nn.Module):
  def __init__(self, in_size=4, depth=4, mlp_depth=4,
               hidden_size=128, attention_size=128,
               relative_size=2, heads=8):
    super().__init__()
    self.good_blocks = nn.ModuleList([
      TransformerBlock(
        relative_size=relative_size,
        hidden_size=hidden_size,
        attention_size=attention_size,
        heads=heads, depth=mlp_depth
      )
      for idx in range(depth)
    ])
    self.bad_blocks = nn.ModuleList([
      TransformerBlock(
        relative_size=relative_size,
        hidden_size=hidden_size,
        attention_size=attention_size,
        heads=heads, depth=mlp_depth
      )
      for idx in range(depth)
    ])

    self.preprocess = nn.Conv1d(in_size, hidden_size, 3, padding=1)
    self.postprocess = nn.Conv1d(hidden_size, 4, 3, padding=1)

  def forward(self, inputs, rna,
              source_bonds, bonds,
              source_stems, stems,
              source_unpaired, unpaired):
    source_structure = DistanceStructure(0, 0, source_bonds.indices, source_bonds.connections, node_count=source_bonds.node_count)
    structure = DistanceStructure(0, 0, bonds.indices, bonds.connections, node_count=bonds.node_count)
    out = scatter.batched(self.preprocess, inputs, rna.indices)
    for good, bad in zip(self.good_blocks, self.bad_blocks):
      out = bad(out, rna, source_structure)
      out = good(out, rna, structure)
    out = scatter.batched(self.postprocess, out, rna.indices)
    return out

class ConvAll(ConvBonds):
  def __init__(self, in_size=5, depth=6, mlp_depth=4,
               hidden_size=128, attention_size=128,
               relative_size=2, heads=8):
    super().__init__(
      in_size=in_size, depth=depth, mlp_depth=mlp_depth,
      hidden_size=hidden_size, attention_size=attention_size,
      relative_size=relative_size, heads=heads
    )
    self.stems = nn.ModuleList([
      TransformerBlock(
        relative_size=relative_size,
        hidden_size=hidden_size,
        attention_size=attention_size,
        heads=heads, depth=mlp_depth
      )
      for idx in range(depth)
    ])
    self.loops = nn.ModuleList([
      TransformerBlock(
        relative_size=relative_size,
        hidden_size=hidden_size,
        attention_size=attention_size,
        heads=heads, depth=mlp_depth
      )
      for idx in range(depth)
    ])

  def forward(self, inputs, rna, bonds, stems, unpaired):
    bonds = DistanceStructure(0, 0, bonds.indices, bonds.connections, node_count=bonds.node_count)
    stems = DistanceStructure(0, 0, stems.indices, stems.connections, node_count=stems.node_count)
    unpaired = DistanceStructure(0, 0, unpaired.indices, unpaired.connections, node_count=unpaired.node_count)

    out = scatter.batched(self.preprocess, inputs, rna.indices)
    for idx, (bond, stem, loop) in enumerate(zip(self.blocks, self.stems, self.loops)):
      out = bond(out, rna, bonds)
      #out = stem(out, rna, stems)
      if idx % 3 == 0:
        out = loop(out, rna, unpaired)
    out = scatter.batched(self.postprocess, out, rna.indices)
    return out
