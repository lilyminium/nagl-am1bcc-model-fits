#!/usr/bin/env bash
#SBATCH -J generate-tetrapeptides
#SBATCH -p free
#SBATCH -t 4:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --account lilyw7
#SBATCH --output slurm-%x.%A.out

. ~/.bashrc

# Use the right conda environment
conda activate openff-nagl-test
conda env export > conda_env.yaml

# Run the commands
python generate-tetrapeptides.py                                \
    --output-smiles-file        "input/tetrapeptides.smi"       \
    --output-fasta-file         "input/tetrapeptides.fasta"     \
    --max-length            4
