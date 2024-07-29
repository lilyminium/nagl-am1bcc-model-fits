#!/usr/bin/env bash
#SBATCH -J label-openeye-am1elf10
#SBATCH --array=1-100%1
#SBATCH -p standard
#SBATCH -t 00:10:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --account dmobley_lab
#SBATCH --output slurm-%x.%A-%a.out

# ===================== conda environment =====================
. ~/.bashrc
conda activate openff-nagl-test


# NAME="all-chembl-small"
# NAME="tetrapeptides"
NAME="small-small-molecules"

echo $NAME


CHARGE_METHOD="am1bccelf10"

python label-openeye-pyarrow.py                                         \
    --input                         "../02_filtered/output/${NAME}"     \
    --output                        "pyarrow/openeye-${CHARGE_METHOD}-${NAME}"  \
    --charge-method                 $CHARGE_METHOD          \
    --n-workers                     600                     \
    --worker-type                   "slurm"                 \
    --batch-size                    200                     \
    --memory                        32                      \
    --walltime                      480                     \
    --queue                         "free"                  \
    --conda-environment             "openff-nagl-test"



python gather-pyarrow-to-datasets.py                \
    --input                         "../02_filtered/output/${NAME}"                 \
    --charges                       "pyarrow/openeye-${CHARGE_METHOD}-${NAME}"      \
    --output                        "output/openeye-${CHARGE_METHOD}-${NAME}"       \
    --column                        "${CHARGE_METHOD}_charges_openeye"

