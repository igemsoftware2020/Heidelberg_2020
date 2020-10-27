# rnasupport

What `torchsupport` does for training classes and `protsupport` does for machine-learning on protein structures, `rnasupport` will do for machine-learning on RNA. That is, `rnasupport` provides model and training primitives for generating training data, as well as training generative models on the task of RNA design. To that end, it provides scripts for classical RNA design, as well as dataset generators and a reinforcement-learning environment for RNA design. To learn more about our RNA-design scripts, take a look at the scripts folder. 

## Setup

Please install `ViennaRNA`, `pytorch` and `torchsupport` in your local python environment. Then run the following in the root directory
of `rnasupport`:

```
python setup.py develop
```
