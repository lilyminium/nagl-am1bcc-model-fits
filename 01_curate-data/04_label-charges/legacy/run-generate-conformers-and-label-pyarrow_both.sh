#!/bin/bash
#
# Set the job name and wall time limit
#BSUB -J label-all-chembl-small-100
#BSUB -W 168:00
#
# Set the output and error output paths.
#BSUB -o  label-all-chembl-small-100-%J.o
#BSUB -e  label-all-chembl-small-100-%J.e
#
#BSUB -M 48

# ===================== conda environment =====================
. ~/.bashrc
micromamba activate openff-nagl-pyarrow

export OE_LICENSE="/home/lilywang/oe_license.txt"


# python generate_conformers_and_label_pyarrow_both.py     \
#     --dataset-file      "../02_filtered/output/enumerated-nci-250k-small-100.smi"   \
#     --output-path       "output/nci-250k-small-100"                                 \
#     --ambertools                \
#     --openeye                   \
#     --n-workers                     300                                         \
#     --worker-type                   "lsf"                   \
#     --batch-size                    20                      \
#     --memory                        24                      \
#     --walltime                      12                      \
#     --queue                         "cpuqueue"              \
#     --conda-environment             "openff-nagl-pyarrow"

python generate_conformers_and_label_pyarrow_both.py     \
    --dataset-file      "../02_filtered/output/enumerated-all-chembl-small-101-199.smi"   \
    --output-path       "output/all-chembl-small-101-199"                                 \
    --ambertools                \
    --openeye                   \
    --n-workers                     300                                         \
    --worker-type                   "lsf"                   \
    --batch-size                    100                     \
    --memory                        24                      \
    --walltime                      12                      \
    --queue                         "cpuqueue"              \
    --conda-environment             "openff-nagl-pyarrow"
