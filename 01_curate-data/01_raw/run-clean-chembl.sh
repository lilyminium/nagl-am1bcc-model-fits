#!/bin/bash
#
# Set the job name and wall time limit
#BSUB -J clean-chembl
#BSUB -W 20:00
#
# Set the output and error output paths.
#BSUB -o  clean-chembl-%J.o
#BSUB -e  clean-chembl-%J.e
#
# Set any cpu options.
#BSUB -M 64


# ===================== conda environment =====================
. ~/.bashrc
micromamba activate openff-nagl-pyarrow

export OE_LICENSE="/home/lilywang/oe_license.txt"


# ===================== run =====================

python clean_smiles.py                                      \
    --n-workers                     300                     \
    --worker-type                   "lsf"                   \
    --batch-size                    1000                    \
    --memory                        24                      \
    --walltime                      12                      \
    --queue                         "cpuqueue"              \
    --conda-environment             "openff-nagl-pyarrow"   \
    --input-file        "input/all-chembl.smi"              \
    --output-file       "output/all-chembl.smi"
