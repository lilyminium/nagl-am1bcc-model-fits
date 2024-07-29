#!/bin/bash
#
# Set the job name and wall time limit
#BSUB -J select-n04
#BSUB -W 168:00
#
# Set the output and error output paths.
#BSUB -o  select-n04-%J.o
#BSUB -e  select-n04-%J.e
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

DATA_INPUT_DIRECTORY="../04_labelled/output"
DATA_OUTPUT_DIRECTORY="output"

mkdir -p $DATA_OUTPUT_DIRECTORY

N_MOL=4
SOURCE_NAME="n04"

# ===================== conda environment =====================
. ~/.bashrc
conda activate openff-nagl-refactor
conda env export > "${LSB_JOBNAME}-environment.yaml"

# ===================== run =====================

export OE_LICENSE="/home/lilywang/oe_license.txt"

SELECTED_OUTPUT_FILE="${DATA_OUTPUT_DIRECTORY}/${SOURCE_NAME}-selected.smi"

python select_diverse_molecules.py                                                              \
        --input-file                "${DATA_INPUT_DIRECTORY}/openeye-enamine-10240.sqlite"      \
        --input-file                "${DATA_INPUT_DIRECTORY}/openeye-enamine-50240.sqlite"      \
        --input-file                "${DATA_INPUT_DIRECTORY}/openeye-chembl.sqlite"             \
        --input-file                "${DATA_INPUT_DIRECTORY}/openeye-zinc.sqlite"               \
        --input-file                "${DATA_INPUT_DIRECTORY}/openeye-nci-250k.sqlite"           \
        --input-file                "${DATA_INPUT_DIRECTORY}/openeye-pdb.sqlite"                \
        --output-file               $SELECTED_OUTPUT_FILE                                       \
        --n-min-molecules           $N_MOL                                                      \
        --element-order             "S"                                                         \
        --element-order             "F"                                                         \
        --element-order             "Cl"                                                        \
        --element-order             "Br"                                                        \
        --element-order             "I"                                                         \
        --element-order             "P"                                                         \
        --element-order             "O"                                                         \
        --element-order             "N"                                                         \
        --element-order             "C"                                                         \
        --element-order             "H"                                                         

