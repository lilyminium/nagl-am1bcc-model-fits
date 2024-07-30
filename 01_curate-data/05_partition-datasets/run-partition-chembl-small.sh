#!/bin/bash
#
# Set the job name and wall time limit
#BSUB -J partition-chembl-small
#BSUB -W 168:00
#
# Set the output and error output paths.
#BSUB -o  partition-chembl-small-%J.o
#BSUB -e  partition-chembl-small-%J.e
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


# ===================== conda environment =====================
. ~/.bashrc
micromamba activate openff-nagl-pyarrow
micromamba env export > "${LSB_JOBNAME}-environment.yaml"

# ===================== run =====================

export OE_LICENSE="/home/lilywang/oe_license.txt"


python partition_diverse_molecules_pyarrow.py                                                           \
        --n-workers                     150                                                                         \
        --worker-type                   "lsf"                                                                       \
        --batch-size                    1000                                                                         \
        --memory                        24                                                                          \
        --walltime                      12                                                                          \
        --queue                         "cpuqueue"                                                                  \
        --conda-environment             "openff-nagl-pyarrow"                                                      \
        --input                     "../04_labelled/output/all-chembl-small-100"                \
        --input                     "../04_labelled/output/all-chembl-small-101-199"            \
        --output-base-name          "output/all-chembl-small"                                   \
        --training-fraction         0.7                                                         

