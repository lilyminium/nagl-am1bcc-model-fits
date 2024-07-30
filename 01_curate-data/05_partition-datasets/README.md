# Partitioning datasets

Datasets are partitioned into training/validation sets here. For `n04_training` and `n04_validation`, this was done following the selection step in `../03_select-diverse`. For the small ChEMBL set, no prior selection step was done.

The partitioning script:

* calculates a Morgan fingerprint as a 2048 bit vector, with radius 3, for each molecule
* uses RDKit's MaxMinPicker to lazily pick diverse molecules by Morgan fingerprint for the training set using the [MaxMin](https://doi.org/10.1002/qsar.200290002) algorithm. This selects the molecule with the highest minimum distance to any molecule already in the picked dataset.
* assigns the remaining molecules to the validation set
* splits the datasets appropriately


## Legacy

A legacy script is provided in `legacy/` that follows the same steps but using the `MoleculeStore` storage format that was used in earlier (v0.1.x, v0.2.x) versions of NAGL.
