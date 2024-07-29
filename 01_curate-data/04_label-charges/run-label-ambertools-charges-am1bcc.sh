#!/usr/bin/env bash
#SBATCH -J label-ambertools
#SBATCH -p standard
#SBATCH -t 72:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --account dmobley_lab
#SBATCH --output slurm-%x.%A.out

# ===================== conda environment =====================
. ~/.bashrc
conda activate openff-nagl-test

export OE_LICENSE=""


NAME="small-small-molecules"

python label-ambertools-pyarrow.py                                      \
    --input                         "../02_filtered/output/${NAME}"     \
    --output                        "pyarrow/ambertools-${NAME}"        \
    --n-workers                     600                     \
    --worker-type                   "slurm"                 \
    --batch-size                    200                     \
    --memory                        32                      \
    --walltime                      480                     \
    --queue                         "free"                  \
    --conda-environment             "openff-nagl-test"



python gather-pyarrow-to-datasets.py                \
    --input                         "../02_filtered/output/${NAME}"     \
    --charges                       "pyarrow/ambertools-${NAME}"        \
    --output                        "output/ambertools-${NAME}"         \
    --column                        "am1bcc_charges_ambertools"
