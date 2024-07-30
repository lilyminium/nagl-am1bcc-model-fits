#!/bin/bash
#
# Set the job name and wall time limit
#BSUB -J prepare-data
#BSUB -W 168:00
#
# Set the output and error output paths.
#BSUB -o  prepare-data-%J.o
#BSUB -e  prepare-data-%J.e
#
#
#BSUB -M 128

# ===================== functions =====================

function check_variable {
    # !1: indirect variable, i.e., expands variable name
    # if $1 value is not set, string is empty and error is printed
    # if it is set, then string == "x"
    test -n "${!1+x}" || echo "$1 not set, "
}

function echovar { echo "${1}='${!1}'" ;}

# ===================== check inputs =====================

source /home/lilywang/.bashrc
micromamba activate openff-nagl-pyarrow

export OE_LICENSE="/home/lilywang/oe_license.txt"

micromamba env export > conda_env.yaml


python prepare_data_pickled.py                                              \
        --n-workers                     500                         \
        --worker-type                   "lsf"                       \
        --batch-size                    1                           \
        --memory                        128                         \
        --walltime                      16                          \
        --queue                         "cpuqueue"                  \
        --conda-environment             "openff-nagl-pyarrow"       \
    --model-config-file         ../01_configs/features/v02.yaml                              \
    --model-config-file         ../01_configs/hyperparameters/reg-am1bcc-0.3.yaml        \
    --data-config-file          ../01_configs/data/am1bcc_small-mol_multi-weighted-objective.yaml          \
    --data-cache-directory      ../01_configs/cached-data                                    \
