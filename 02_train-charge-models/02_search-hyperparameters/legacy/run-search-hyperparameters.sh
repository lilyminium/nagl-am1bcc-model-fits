#!/bin/bash
#
# Set the job name and wall time limit
#BSUB -J search-hyperparameters
#BSUB -W 168:00
#
# Set the output and error output paths.
#BSUB -o  search-hyperparameters-%J.o
#BSUB -e  search-hyperparameters-%J.e
#
# Set any gpu options.
#BSUB -q gpuqueue
#BSUB -gpu num=1:j_exclusive=yes:mode=shared:mps=no:
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

variable_error=$(
    check_variable "CHARGE_METHOD"
    check_variable "FEATURE_NAME"
    check_variable "DATASET_NAME"
    check_variable "HYPERPARAMETER_NAME"
)
echo $variable_error
test -n "${variable_error}" && exit 1


# ===================== environment =====================

CONFIG_DIRECTORY="../01_configs"
CACHE_DIRECTORY="../00_cached-data"
OUTPUT_DIRECTORY="./output"
mkdir -p $OUTPUT_DIRECTORY


# ===================== conda environment =====================
. ~/.bashrc
conda activate openff-nagl

conda env export > "${OUTPUT_DIRECTORY}/${LSB_JOBNAME}-environment.yaml"


FEATURE_CONFIG="${CONFIG_DIRECTORY}/features/${FEATURE_NAME}.yaml"
DATA_CONFIG="${CONFIG_DIRECTORY}/data/${DATASET_NAME}.yaml"
OUTPUT_MODEL_CONFIG="${CONFIG_DIRECTORY}/hyperparameters/${HYPERPARAMETER_NAME}.yaml"

echovar "FEATURE_CONFIG"
echovar "DATA_CONFIG"
echovar "OUTPUT_MODEL_CONFIG"

echovar "CACHE_DIRECTORY"
echovar "CHARGE_METHOD"


python search_hyperparameters.py                                            \
    --model-config-file         $FEATURE_CONFIG                             \
    --model-config-file         $DATA_CONFIG                                \
    --output-directory          $OUTPUT_DIRECTORY                           \
    --n-epochs                  250                                         \
    --n-total-trials            400                                         \
    --n-gpus                    1                                           \
    --partial-charge-method     $CHARGE_METHOD                              \
    --data-cache-directory      $CACHE_DIRECTORY                            \
    --output-config-file        $OUTPUT_MODEL_CONFIG                        