#!/bin/bash
#

# wget http://ligand-expo.rcsb.org/dictionaries/Components-smiles-stereo-oe.smi
# mv Components-smiles-stereo-oe.smi output/pdb.smi

# curl "https://zinc.docking.org/substances/subsets/fda.smi?count=all" > input/fda.smi
cat input/fda.smi | awk '{print $1}' > output/fda.smi