#!/usr/bin/env bash
#SBATCH -J select-n04
#SBATCH -p free
#SBATCH -t 4:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --account lilyw7
#SBATCH --output slurm-%x.%A.out

# ===================== conda environment =====================
. ~/.bashrc
conda activate openff-nagl-test

# ===================== run =====================


DATA_DIRECTORY="../02_filtered/output"


python select-from-fingerprints.py                              \
    --input             "${DATA_DIRECTORY}/enamine-10240"       \
    --input             "${DATA_DIRECTORY}/enamine-50240"       \
    --input             "${DATA_DIRECTORY}/chembl"              \
    --input             "${DATA_DIRECTORY}/zinc"                \
    --input             "${DATA_DIRECTORY}/nci-250k"            \
    --input             "${DATA_DIRECTORY}/pdb"                 \
    --n-min-molecules   4                                       \
    --output            "output/n04-all"
