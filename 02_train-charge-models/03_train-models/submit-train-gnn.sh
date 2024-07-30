#!/bin/bash


# train v0.1.0-rc1

bsub -env CHARGE_METHOD="am1bcc",FEATURE_NAME="v02",OBJECTIVE="multi",DATASET="n04-peptides",REPLICATE="rep-01" < run-train-reg.sh

# train v0.1.0-rc2
bsub -env CHARGE_METHOD="am1bcc",FEATURE_NAME="v02",OBJECTIVE="multi",DATASET="small-mol",REPLICATE="rep-01" < run-train-reg-from-start.sh

