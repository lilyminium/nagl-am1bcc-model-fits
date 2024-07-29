#!/bin/bash
#
# Set the job name and wall time limit
#BSUB -J label-openeye-enamine-10240
#BSUB -W 168:00
#
# Set the output and error output paths.
#BSUB -o  label-openeye-enamine-10240-%J.o
#BSUB -e  label-openeye-enamine-10240-%J.e
#
#BSUB -M 48

# ===================== conda environment =====================
. ~/.bashrc
conda activate openff-nagl-refactor

# ===================== run =====================
SMILES_INPUT_DIRECTORY="../02_filtered/output"

export OE_LICENSE="/home/lilywang/oe_license.txt"

python generate_conformers_and_label.py                                                                         \
    --smiles-file                   "${SMILES_INPUT_DIRECTORY}/enumerated-mapped-enamine-10240.smi"                            \
    --output-file                   "output/openeye-enamine-10240.sqlite"                                               \
    --partial-charge-method         "am1bcc"                                                                    \
    --partial-charge-method         "am1"                                                                       \
    --n-workers                     200                                                                         \
    --worker-type                   "lsf"                                                                       \
    --batch-size                    5                                                                           \
    --memory                        4                                                                          \
    --walltime                      24                                                                          \
    --queue                         "cpuqueue"                                                                  \
    --conda-environment             "openff-nagl-refactor"                                                   
