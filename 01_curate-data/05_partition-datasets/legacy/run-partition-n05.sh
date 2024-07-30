#!/bin/bash
#
# Set the job name and wall time limit
#BSUB -J partition-n05
#BSUB -W 168:00
#
# Set the output and error output paths.
#BSUB -o  partition-n05-%J.o
#BSUB -e  partition-n05-%J.e
#
# Set any cpu options.
#BSUB -n 1 -R "span[ptile=1]"
#BSUB -M 128

# ===================== functions =====================

function check_variable {
    # !1: indirect variable, i.e., expands variable name
    # if $1 value is not set, string is empty and error is printed
    # if it is set, then string == "x"
    test -n "${!1+x}" || echo "$1 not set, "
}

function echovar { echo "${1}='${!1}'" ;}


# ===================== environment =====================

N_MOL=5
SOURCE_NAME="n05"

# ===================== conda environment =====================
. ~/.bashrc
conda activate openff-nagl-refactor
conda env export > "${LSB_JOBNAME}-environment.yaml"

# ===================== run =====================

export OE_LICENSE="/home/lilywang/oe_license.txt"

SELECTED_OUTPUT_FILE="../03_select-diverse/output/${SOURCE_NAME}-selected.smi"
OUTPUT_BASE_NAME="input/${SOURCE_NAME}"

python partition_diverse_molecules.py                                                           \
        --n-workers                     150                                                                         \
        --worker-type                   "lsf"                                                                       \
        --batch-size                    1000                                                                         \
        --memory                        24                                                                          \
        --walltime                      12                                                                          \
        --queue                         "cpuqueue"                                                                  \
        --conda-environment             "openff-nagl-refactor"                                                      \
        --input-file                $SELECTED_OUTPUT_FILE                                       \
        --output-base-name          $OUTPUT_BASE_NAME                                           \
        --training-fraction         0.8                                                         

