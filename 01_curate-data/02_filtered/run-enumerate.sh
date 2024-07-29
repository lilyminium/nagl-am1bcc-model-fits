#!/usr/bin/env bash
#SBATCH -J enumerate
#SBATCH -p free
#SBATCH -t 4:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --account lilyw7
#SBATCH --output slurm-%x.%A.out


# ===================== conda environment =====================
. ~/.bashrc
conda activate openff-nagl-test

# ===================== run =====================


python enumerate.py                                                     \
    --input-file            "../01_raw/input/tetrapeptides.smi"         \
    --output-file           "enumerated/enumerated-tetrapeptides.smi"   \
    --no-tautomers                                                      \
    --protomers                                                         \
    --max-protomers         2                                           \
    --n-processes           4

