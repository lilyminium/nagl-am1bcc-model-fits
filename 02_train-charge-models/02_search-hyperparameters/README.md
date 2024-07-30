# Hyperparameter searching

Early on in the project we did a quick search over hyperparameters using Ray Tune. Code for this only exists for the 0.2.x version of NAGL, in the `legacy` directory.

Hyperparameters were searched on a smaller subset of the data, using just the OpenFF Industry Benchmark set as both the training and validation sets. As this is before multiple targets were supported, the training target is only fitting to AM1-BCC charges.

The resulting "best" hyperparameters for each configuration searched are collated in the `../01_configs/legacy/hyperparameters` directory.
We searched the best hyperparameters for reproducing both "am1" and "am1bcc" charges, across a number of features. These are listed below, along with the new name of the feature set in the current `../01_configs/features` directory.


The following hyperparameters were searched:

* number of graph convolution layers: 4, 5, 6
* number of hidden features per graph convolution layer: 128, 256, 512
* graph convolution architecture: SAGEConv, GINConv
* number of dense layers: 1, 2
* number of hidden features per dense layer: 64, 128, 256
* activation function: ReLU, Sigmoid
* learning rate: 5e-2, 1e-3, 5e-3


The most common features in each category were then collated in the "best-am1-common.yaml" and "best-am1bcc-common.yaml" files. This was done manually in the `../01_curate-data/legacy/collated-*-hyperparameters.csv` files.

