# Filtering and enumerating tautomers

The scripts in this directory handle filtering molecules from datasets, as well as enumerating protomers or tautormers.
Typically a dataset is first filtered before being enumerated.

This directory also contains a script for converting SMILES into [Hugging Face datasets](https://huggingface.co/docs/datasets/en/index).
However, this step is very optional. While some scripts later on (e.g. in `../03_select-diverse`) expect a Dataset input, it is trivial to re-write them to simply use SMILES.

## Filtering

This is done using the `filter_smiles.py`, or in versions up to 0.2.x NAGL was done using the `filter_molecules`
function in `openff.nagl._app`. This allows filtering for:

* minimum and maximum mass
* number of rotatable bonds
* allowed elements
* excluded SMILES (e.g. those in the test dataset)


## Enumerating protomers and/or tautomers

This script gives the option of enumerating protomers or tautomers, using the OpenFF or OpenEye toolkits.
Typically for each dataset we enumerated protomers (up to a max of 2), but not tautomers.
