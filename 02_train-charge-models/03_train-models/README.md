# Training files

The files here take on two forms:

1. Pickling training data

The script `prepare_data_pickled.py` prepares the training data and pickles it in a parallelized way for efficiency. This is an optional step that is not necessary for smaller datasets, which should get processed by NAGL as it sets up the `TrainingGNNModel`. Data is stored in the `01_configs/cached-data` directory.

2. Training the neural network

This is accomplished with `train_gnn.py`, where `run-train-reg.sh` shows an example use of the Python script being invoked. `submit-train-gnn.sh` shows submitting `run-train-reg.sh` with particular environment variables.


### Loading the trained networks

Outputs are stored in the `output` directory. `best_model.pt` is a `GNNModel` that can be loaded with:

```python
from openff.nagl import GNNModel

model = GNNModel.load("best_model.pt", eval_mode=True)
```

If using a version of NAGL that requires version numbers, use the file instead named `output_*.pt` -- these were created with `export-to-versioned-model.py`. The export script simply adds a version number to models that were trained before that property existed.


### Model history

The NAGL 0.1.0-rc1 model is from `output/am1bcc_v02_n04-peptides_multi_reg/rep-01/output_am1bcc_v02_n04-peptides_multi_reg_rep-01.pt`.
This was trained using `run-train-reg.sh` with the following config files:

* features/v02.yaml
* data/am1bcc_n04-peptides_multi-objective.yaml
* hyperparameters/reg-am1bcc-0.3.yaml

The NAGL 0.1.0-rc2 model is from `output/am1bcc_v02_small-mol_multi_reg/rep-01/best_model.pt`. This was trained using NAGL 0.1.0-rc1 as a starting point and training further to a small molecule dataset.

This was done using `run-train-reg-from-start.sh` with the following config files:

* features/v02.yaml
* data/am1bcc_small-mol_multi-objective.yaml
* hyperparameters/reg-am1bcc-0.3.yaml


