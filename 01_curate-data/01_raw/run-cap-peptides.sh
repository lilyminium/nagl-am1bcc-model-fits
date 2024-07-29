#!/bin/bash
#SBATCH -J cap-peptides
#SBATCH --array=0-3
#SBATCH -p free
#SBATCH -t 4:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --account lilyw7
#SBATCH --output slurm-%x.%A-%a.out



# ===================== conda environment =====================
. ~/.bashrc
conda activate openff-nagl-test

# create array of peptide chain lengths
# N_PEPTIDE_CHAINS=(5 6 7 8 9 10 12 15 18 21 30)
N_PEPTIDE_CHAINS=(1 2 3 4)

# index with slurm job index
i=${N_PEPTIDE_CHAINS[$SLURM_ARRAY_TASK_ID]}

echo $i

python cap-peptides.py             \
    --input-file            "input/peptide-chains-${i}.fasta"              \
    --output-smiles-file    "input/peptide-chains-${i}.smi"                \
    --output-mapped-smiles-file     "input/peptide-chains-${i}-mapped.smi" \

