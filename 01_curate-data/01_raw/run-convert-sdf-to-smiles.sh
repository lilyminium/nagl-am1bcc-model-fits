#!/bin/bash
#SBATCH -J convert-sdf
#SBATCH -p free
#SBATCH -t 4:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --account lilyw7
#SBATCH --output slurm-%x.%A.out



# ===================== conda environment =====================
. ~/.bashrc
conda activate openff-nagl-test
conda env export > "${LSB_JOBNAME}-environment.yaml"

python -c "import openff.nagl ; print(openff.nagl.__version__)"


# ===================== run =====================

mkdir -p output

# export OE_LICENSE="/home/lilywang/oe_license.txt"

python convert_sdf_to_smiles.py                                                     \
    --sdf-file              "input/enamine-10240.sdf.gz"                            \
    --output-file           "output/enamine-10240.smi"                              
