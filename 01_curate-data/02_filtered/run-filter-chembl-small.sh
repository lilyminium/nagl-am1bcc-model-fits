#!/bin/bash
#
# Set the job name and wall time limit
#BSUB -J filter-all-chembl-small
#BSUB -W 20:00
#
# Set the output and error output paths.
#BSUB -o  filter-all-chembl-small-%J.o
#BSUB -e  filter-all-chembl-small-%J.e
#
# Set any cpu options.
#BSUB -M 64
#BSUB -n 8



# ===================== conda environment =====================
. ~/.bashrc
micromamba activate openff-nagl-refactor
micromamba env export > "${LSB_JOBNAME}-environment.yaml"


# ===================== run =====================

mkdir -p output input

export OE_LICENSE="/home/lilywang/oe_license.txt"

SMILES_DIRECTORY="../01_raw/output"


for dataset in "all-chembl" ; do

    python filter_smiles.py                                                 \
        --input-file            "${SMILES_DIRECTORY}/${dataset}.smi"        \
        --output-file           "input/filtered-${dataset}-small.smi"       \
        --exclude-smiles        "output/enumerated-benchmark.smi"           \
        --exclude-smiles        "output/enumerated-spice.smi"               \
        --min-mass              1                                           \
        --max-mass              199                                         \
        --n-rotatable-bonds     10                                          \
        --only-retain-largest                                               \
        --n-processes           8

done

