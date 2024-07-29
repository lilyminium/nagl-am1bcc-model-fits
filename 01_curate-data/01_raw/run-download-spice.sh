#!/bin/bash

. ~/.bashrc
conda activate openff-nagl

wget https://github.com/choderalab/charge-datasets/releases/download/1.1.0/spice.oeb

python extract_smiles_from_oeb.py                       \
    --input-file        spice.oeb                       \
    --output-file       output/spice.smi