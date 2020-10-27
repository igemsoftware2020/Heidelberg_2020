from RNA import fold

import torch
import torch.nn as nn
import torch.nn.functional as func

class ConnectedConv(nn.Module):
  def __init__(self, in_size, out_size, kernel_size=3, position_size=10,
               causal=False, dilation=1):
    super().__init__()
    self.causal = causal
    self.position_size = position_size
    self.kernel_size = kernel_size
    self.dilation = dilation
    self.linear = nn.Conv1d(
      (2 * in_size + position_size) * kernel_size,
      out_size, 1
    )

  def position_encoding(self, connections):
    positions = torch.arange(connections.size(1))[None].float()
    positions = positions.to(connections.device)
    return torch.cat([
      torch.sin(2 ** idx * (positions - connections.float()) / 1000)[:, None]
      for idx in range(self.position_size)
    ], dim=1)

  def forward(self, inputs, connections, mask):
    ind = torch.arange(inputs.size(0), dtype=torch.long, device=inputs.device)
    connection_values = inputs[ind[:, None], :, connections].transpose(1, 2)
    position = self.position_encoding(connections)
    out = torch.cat((inputs, connection_values, position), dim=1)
    if self.causal:
      padding = (self.kernel_size * self.dilation, 0)
    else:
      padding = (
        self.kernel_size * self.dilation // 2,
        self.kernel_size * self.dilation // 2
      )
    out = func.pad(out, padding)
    out = out.unfold(-1, self.kernel_size, 1)
    if self.causal:
      out = out[:, :, :-1]
    out = out.permute(0, 1, 3, 2).reshape(inputs.size(0), -1, inputs.size(2))
    out = self.linear(out)
    out[~mask[:, None, :].expand_as(out)] = 0
    return out

class ConnectedConvBlock(nn.Module):
  def __init__(self, size, kernel_size=5, position_size=10, causal=True):
    super().__init__()
    self.block = ConnectedConv(
      size, size, kernel_size=kernel_size,
      position_size=position_size, causal=causal
    )
    self.mix = nn.Conv1d(size, size, 1)
    self.bn = nn.InstanceNorm1d(size)

  def forward(self, inputs, connections, mask):
    out = self.block(self.bn(inputs), connections, mask)
    out = func.relu(out)
    out = self.mix(out)
    out = func.relu(out + inputs)
    return out

class ConnectedConvEncoder(nn.Module):
  def __init__(self, size, kernel_size=5, position_size=10, depth=4):
    super().__init__()
    self.position_size = position_size
    self.predict = nn.Conv1d(position_size, size, 1)
    self.blocks = nn.ModuleList([
      ConnectedConvBlock(
        position_size,
        kernel_size=kernel_size,
        position_size=position_size,
        causal=False
      )
      for idx in range(depth)
    ])

  def position_encoding(self, connections):
    positions = torch.arange(connections.size(1))[None].float()
    positions = positions.to(connections.device)
    result = torch.cat([
      torch.sin(2 ** idx * positions / 1000)[:, None]
      for idx in range(self.position_size)
    ], dim=1)
    return result.repeat_interleave(connections.size(0), dim=0)

  def forward(self, connections, mask):
    out = self.position_encoding(connections)
    for block in self.blocks:
      out = block(out, connections, mask)
    out = self.predict(out)
    out[~mask[:, None].expand_as(out)] = 0
    return out

class ConnectedConvDecoder(nn.Module):
  def __init__(self, size, kernel_size=5, position_size=10, depth=4):
    super().__init__()
    self.position_size = position_size
    self.encode = ConnectedConv(
      4, size, kernel_size=1,
      position_size=position_size,
      causal=True
    )
    self.blocks = nn.ModuleList([
      ConnectedConvBlock(
        2 * size,
        kernel_size=kernel_size,
        position_size=position_size,
        causal=True
      )
      for idx in range(depth)
    ])
    self.predict = nn.Conv1d(2 * size, 4, 1)

  def forward(self, inputs, code, connections, mask):
    out = torch.cat((self.encode(inputs, connections, mask), code), dim=1)
    for block in self.blocks:
      out = block(out, connections, mask)
    return self.predict(out)

class ConnectedConvTranscoder(nn.Module):
  def __init__(self, size, kernel_size=5, position_size=10, depth=4):
    super().__init__()
    self.position_size = position_size
    self.encode = nn.Conv1d(4, size, 1)
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
    self.predict = nn.Conv1d(size, 4, 1)
    self.position = nn.Conv1d(size, 1, 1)

  def forward(self, sequence, current, target, mask):
    out = self.encode(sequence)
    for cblock, tblock in zip(self.current_blocks, self.target_blocks):
      out = cblock(out, current, mask)
      out = tblock(out, target, mask)
    out = self.predict(out)
    out[~mask[:, None].expand_as(out)] = 0
    return out

class ConnectedConvDiffusion(nn.Module):
  def __init__(self, size, kernel_size=5, position_size=10, depth=4):
    super().__init__()
    self.position_size = position_size
    self.encode = nn.Conv1d(4 + position_size, size, 1)
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
    self.predict = nn.Conv1d(size, 4, 1)
    self.position = nn.Conv1d(size, 1, 1)

  def forward(self, sequence, step, current, target, mask):
    step = step[:, :, None].repeat_interleave(sequence.size(2), dim=2)
    out = self.encode(torch.cat((sequence, step), dim=1))
    for cblock, tblock in zip(self.current_blocks, self.target_blocks):
      out = cblock(out, current, mask)
      out = tblock(out, target, mask)
    out = self.predict(out)
    out[~mask[:, None].expand_as(out)] = 0
    return out

class ConnectedConvModel(nn.Module):
  def __init__(self, size, kernel_size=5, position_size=10, depth=4):
    super().__init__()
    self.encoder = ConnectedConvEncoder(
      size, kernel_size=kernel_size,
      position_size=position_size,
      depth=depth
    )
    self.decoder = ConnectedConvDecoder(
      size, kernel_size=kernel_size,
      position_size=position_size,
      depth=depth
    )

  def forward(self, inputs, bonds, causal_bonds, mask):
    code = self.encoder(bonds, mask)
    prediction = self.decoder(inputs, code, causal_bonds, mask)
    return prediction

class ConnectedConvRL(ConnectedConvModel):
  def sample(self, bonds, causal_bonds, mask):
    code = self.encoder(bonds, mask)
    out = torch.zeros(mask.size(0), 4, mask.size(1)).to(mask.device)
    predictions = []
    inds = torch.arange(mask.size(0), dtype=torch.long, device=mask.device)
    for idx in range(mask.size(1)):
      pred = self.decoder(out, code, causal_bonds, mask)
      pred = func.log_softmax(pred, dim=1)
      distribution = torch.distributions.OneHotCategorical(
        logits=pred.transpose(1, 2)
      )
      samples = distribution.sample().transpose(1, 2)
      out[:, :, idx] = samples[:, :, idx]
      pos = samples[:, :, idx].argmax(dim=1)
      predictions.append(pred[inds, pos, idx].unsqueeze(1))
    predictions = torch.cat(predictions, dim=1)
    seqs = out.argmax(dim=1)
    return predictions, seqs

  def discount(self, reward, mask):
    count = sum(mask.long())
    result = torch.zeros(mask.size(0), device=mask.device)
    for idx in range(count):
      result[count - idx] = reward * 0.9 ** idx
    return result

  def reward(self, sequences, target, mask):
    result = []
    for seq, t, m in zip(sequences, target, mask):
      struc, _ = fold("".join(["GAUC"[x] for x, mm in zip(seq, m) if mm]))
      t = "".join([".()"[x] for x, mm in zip(t, m) if mm])
      reward = int(struc == t)
      result.append(self.discount(reward, m).unsqueeze(0))
    return torch.cat(result, dim=0)

  def rl_sample(self, structure, bonds, causal_bonds, mask):
    predictions, seqs = self.sample(bonds, causal_bonds, mask)
    rewards = self.reward(seqs, structure, mask)
    return predictions, rewards

  def forward(self, inputs, structure, bonds, causal_bonds, mask):
    teacher_forced = super().forward(inputs, bonds, causal_bonds, mask)
    reinforced = self.rl_sample(structure, bonds, causal_bonds, mask)
    return teacher_forced, reinforced

class ConnectedConvScheduled(ConnectedConvModel):
  def __init__(self, size, steps=10, mixing=0.5, **kwargs):
    super().__init__(size, **kwargs)
    self.steps = steps
    self.mixing = mixing

  def forward(self, inputs, bonds, causal_bonds, mask):
    with torch.no_grad():
      for idx in range(self.steps - 1):
        predictions = super().forward(inputs, bonds, causal_bonds, mask)
        distribution = torch.distributions.OneHotCategorical(
          logits=predictions.transpose(1, 2)
        )
        samples = distribution.sample().transpose(1, 2)
        mask = torch.rand(inputs.size(0), inputs.size(2)) < self.mixing
        inputs.transpose(1, 2)[mask] = samples.transpose(1, 2)[mask]
    return super().forward(inputs, bonds, causal_bonds, mask)
