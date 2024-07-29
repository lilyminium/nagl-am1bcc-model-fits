#!/usr/bin/env bash
#SBATCH -J generate-small-small-molecules
#SBATCH -p free
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128gb
#SBATCH --account lilyw7
#SBATCH --output slurm-%x.%A.out

. ~/.bashrc

# Use the right conda environment
conda activate openff-nagl-test
conda env export > conda_env.yaml

# Run the commands
python generate-small-small-molecules.py > small-small-molecules.log 2>&1
