#!/bin/bash
#
# Set the job name and wall time limit
#BSUB -J train-gnn-500-reg
#BSUB -W 168:00
#
# Set the output and error output paths.
#BSUB -o  train-gnn-500-reg-%J.o
#BSUB -e  train-gnn-500-reg-%J.e
#
# Set any gpu options.
#BSUB -q gpuqueue
#BSUB -gpu num=1:j_exclusive=yes:mode=shared:mps=no:
#
#BSUB -M 496

# ===================== functions =====================

function check_variable {
    # !1: indirect variable, i.e., expands variable name
    # if $1 value is not set, string is empty and error is printed
    # if it is set, then string == "x"
    test -n "${!1+x}" || echo "$1 not set, "
}

function echovar { echo "${1}='${!1}'" ;}

# ===================== check inputs =====================

variable_error=$(
    check_variable "CHARGE_METHOD"
    check_variable "FEATURE_NAME"
    check_variable "OBJECTIVE"
    check_variable "DATASET"
    check_variable "REPLICATE"
)
echo $variable_error
test -n "${variable_error}" && exit 1

# ===================== conda environment =====================

source /home/lilywang/.bashrc
micromamba activate openff-nagl-pyarrow

export OE_LICENSE="/home/lilywang/oe_license.txt"

# CHARGE_METHOD="am1bcc"
# FEATURE_NAME="v02"
# OBJECTIVE="multi"
# DATASET="n04"
# REPLICATE="rep-01"

OUTPUT_SUBDIRECTORY="${CHARGE_METHOD}_${FEATURE_NAME}_${DATASET}_${OBJECTIVE}"
OUTPUT_DIRECTORY="output/${OUTPUT_SUBDIRECTORY}_reg/${REPLICATE}"

echo "OUTPUT_DIRECTORY=${OUTPUT_DIRECTORY}"
mkdir -p $OUTPUT_DIRECTORY

python train_gnn.py                                                 \
    --model-config-file         "../01_configs/features/${FEATURE_NAME}.yaml"                              \
    --model-config-file         "../01_configs/hyperparameters/reg-${CHARGE_METHOD}-0.3.yaml"        \
    --data-config-file          "../01_configs/data/${CHARGE_METHOD}_${DATASET}_${OBJECTIVE}-objective.yaml"          \
    --output-directory          $OUTPUT_DIRECTORY                           \
    --n-epochs                  500                                         \
    --n-gpus                    1                                           \
    --data-cache-directory      ../01_configs/cached-data                            
