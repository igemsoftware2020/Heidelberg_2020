#****************************************************
# Language Model for Protein Generation
# ***************************************************

# Source of the code: https://github.com/salesforce/ctrl/blob/master/pytorch_transformer.py


import torch
import torch.nn as nn
import numpy as np

def freeze_params(encoder):
  for param in encoder.parameters():
    param.requires_grad = False
    
def unfreeze_layer(layer):
  for param in layer.parameters():
    param.requires_grad = True
    
def subsequent_mask(size):
    "Mask out subsequent positions."
    attn_shape = (1, size, size)
    subsequent_mask = np.triu(np.ones(attn_shape), k=1).astype('uint8')
    return ~(torch.from_numpy(subsequent_mask) == 0)
  
# Positional Encoding

def angle_defn(pos, i, d_model_size):
  angle_rates = 1 / np.power(10000, (2 * (i//2)) / np.float32(d_model_size))
  return pos * angle_rates

def positional_encoding(position, d_model_size):
  angle_rads = angle_defn(np.arange(position)[:, np.newaxis], np.arange(d_model_size)[np.newaxis, :], d_model_size)
  
  sines = np.sin(angle_rads[:, 0::2])
  cosines = np.cos(angle_rads[:, 1::2])
  
  pos_encoding = torch.tensor(np.concatenate([sines, cosines], axis=-1)[np.newaxis, ...], dtype = torch.float)
  return pos_encoding

# Calculating the scaled dot product attention
def scaled_dot_product_attention(q, k, v, mask):
  """
  Computes the scaled dot product attention. It returns an attention matrix.

  Args:
    q: queries
    k: keys
    v: values
    mask: mask indicating elements which should not included for the computation.
  """

  matmul_qk = torch.matmul(q, k.permute(0,1,3,2))

  dk = k.shape[-1]
  
  scaled_attention_logits = matmul_qk / np.sqrt(dk)

  if mask is not None:
  
    mask = mask[:, None].repeat_interleave(matmul_qk.size(1), dim=1)
    scaled_attention_logits += (mask * -1e9)

  
  attention_weights = torch.softmax(scaled_attention_logits, dim=-1)

  output = torch.matmul(attention_weights, v)
  
  return output

# Multihead attention
class MultiHeadAttention(torch.nn.Module):
  def __init__(self, d_model_size, num_heads):
    super(MultiHeadAttention, self).__init__()

    self.num_heads = num_heads
    self.d_model_size = d_model_size

    self.depth = int(d_model_size / self.num_heads)
    self.Wq = torch.nn.Linear(d_model_size, d_model_size)
    self.Wk = torch.nn.Linear(d_model_size, d_model_size)
    self.Wv = torch.nn.Linear(d_model_size, d_model_size)
    
    self.dense = torch.nn.Linear(d_model_size, d_model_size)
        
  def split_into_heads(self, x, batch_size):
    x = x.reshape(batch_size, -1, self.num_heads, self.depth)
    return x.permute([0, 2, 1, 3])
    
  def forward(self, v, k, q, mask):
    batch_size = q.shape[0]
    q = self.Wq(q)
    k = self.Wk(k)
    v = self.Wv(v)
    
    q = self.split_into_heads(q, batch_size)
    k = self.split_into_heads(k, batch_size)
    v = self.split_into_heads(v, batch_size)
    
    scaled_attention = scaled_dot_product_attention(q, k, v, mask).permute([0, 2, 1, 3])
    original_size_attention = scaled_attention.reshape(batch_size, -1, self.d_model_size)
    output = self.dense(original_size_attention)
        
    return output


# Fully connected feed forward networkt
def point_wise_feed_forward_network(d_model_size, dff):
  return torch.nn.Sequential(torch.nn.Linear(d_model_size, dff), torch.nn.ReLU(), torch.nn.Linear(dff, d_model_size))


# One Encoder Layer
class EncoderLayer(torch.nn.Module):
  def __init__(self, d_model_size, num_heads, dff, rate=0.1):
    super(EncoderLayer, self).__init__()

    # Initilize the components of the encoder layer
    self.multi_head_attention = MultiHeadAttention(d_model_size, num_heads)
    self.activation = torch.nn.ReLU()
    self.linear1 = torch.nn.Linear(d_model_size, dff)
    self.linear2 = torch.nn.Linear(dff, d_model_size)
    self.rate = torch.nn.Dropout(rate)
    self.dropout1 = torch.nn.Dropout(rate)
    self.dropout2 = torch.nn.Dropout(rate)
    
    # residual weight/ /alpha_i: is initilized with zero
    self.resweight = torch.nn.Parameter(torch.tensor(0.0))

    
  def forward(self, x, mask):

    seq = x
    seq  = self.multi_head_attention(seq, seq, seq, mask)

    seq = seq[0]
    seq = seq * self.resweight
    # xi+1 = xi + resweight*sublayer(xi)
    out1 = x + self.dropout1(seq)

    seq = out1
    seq = self.linear2(self.rate(self.activation(self.linear1(seq))))
    seq = seq * self.resweight
    # xi+1 = xi + resweight*sublayer(xi)
    out2 = out1 + self.dropout2(seq) 
    
    return out2


# Encoder

class Encoder(torch.nn.Module):
  def __init__(self, num_layers=9, d_model_size=512, num_heads=8, dff=512, input_vocab_size=21,
               rate=0.1, **kwargs):
    super(Encoder, self).__init__()

    self.d_model_size = d_model_size
    self.num_layers = num_layers
    
    self.pos_encoding = positional_encoding(self.d_model_size, self.d_model_size)

    for i in range(num_layers):
      setattr(self, "layer%i" % i, EncoderLayer(d_model_size, num_heads, dff, rate))

    self.verkleinern = nn.Conv1d(2*d_model_size, d_model_size, 1)
    self.embedding = nn.Embedding(input_vocab_size, d_model_size)
    self.dropout = torch.nn.Dropout(rate)
    self.predict = nn.Conv1d(d_model_size, input_vocab_size, 1)

  def forward(self, seq, anno_vec, mask):

    seq = self.embedding(seq)

    anno_vec = anno_vec
    mask = mask
    mask = torch.repeat_interleave(mask.unsqueeze(1), 512, dim = 1).bool()
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    mask_fw = subsequent_mask(512)
    mask_fw = mask_fw.to(device)
    mask = mask + mask_fw
    mask = mask.float()

    seq *= np.sqrt(self.d_model_size)
    seq += self.pos_encoding.to(seq.device)
    seq = self.dropout(seq)


    anno_vec = torch.repeat_interleave(anno_vec.unsqueeze(-1), 512, dim = 2)
    

    seq = torch.cat((seq.float(), anno_vec.float()), dim = 1)
    seq = self.verkleinern(seq)
    
    for i in range(self.num_layers):
      seq = getattr(self, "layer%i" % i)(seq, mask)
    out = seq
    return self.predict(out)
