#!/usr/bin/env bash
#SBATCH -J generate-diastereomers
#SBATCH -p standard
#SBATCH -t 48:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --account dmobley_lab
#SBATCH --output slurm-%x.%A.out

# ===================== conda environment =====================
. ~/.bashrc
conda activate openff-nagl-test


python generate-diastereomers.py                                        \
        --n-workers                     500                     \
        --worker-type                   "slurm"                 \
        --batch-size                    10                      \
        --memory                        48                      \
        --walltime                      480                     \
        --queue                         "free"                  \
        --conda-environment             "openff-nagl-test"      \
