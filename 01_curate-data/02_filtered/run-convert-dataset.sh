#!/usr/bin/env bash
#SBATCH -J convert-dataset
#SBATCH -p free
#SBATCH -t 4:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --account lilyw7
#SBATCH --output slurm-%x.%A.out

# ===================== conda environment =====================
. ~/.bashrc
conda activate openff-nagl-test

# ===================== run =====================

# python convert-to-dataset.py                                            \
#     --input             "enumerated/enumerated-all-chembl-small.smi"    \
#     --dataset           "CHEMBL, small"                                 \
#     --output            "output/all-chembl-small"                       \
#     --n-processes       4

python convert-to-dataset.py                                        \
    --input             "enumerated/enumerated-tetrapeptides.smi"   \
    --dataset           "Tetrapeptides"                             \
    --output            "output/tetrapeptides"                      \
    --n-processes       4