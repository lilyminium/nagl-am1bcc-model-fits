#!/bin/bash

# bsub -env CHARGE_METHOD="am1",DATASET_NAME="partitioned-nb-h-n10",FEATURE_NAME="v1-simple",HYPERPARAMETER_NAME="best-am1-v1" < run-search-hyperparameters.sh
# 12073186
# killed from lack of time / slowness

# bsub -env CHARGE_METHOD="am1bcc",DATASET_NAME="partitioned-nb-h-n10",FEATURE_NAME="v1-simple",HYPERPARAMETER_NAME="best-am1bcc-v1" < run-search-hyperparameters.sh
# 11895985
# killed from lack of time / slowness

# bsub -env CHARGE_METHOD="am1",DATASET_NAME="partitioned-nb-h-n10",FEATURE_NAME="v3-espaloma",HYPERPARAMETER_NAME="best-am1-v3" < run-search-hyperparameters.sh
# 11895984
# bsub -env CHARGE_METHOD="am1bcc",DATASET_NAME="partitioned-nb-h-n10",FEATURE_NAME="v3-espaloma",HYPERPARAMETER_NAME="best-am1bcc-v3" < run-search-hyperparameters.sh
# 11895986

# bsub -env CHARGE_METHOD="am1",DATASET_NAME="openff-benchmark",FEATURE_NAME="v1-simple",HYPERPARAMETER_NAME="best-am1-v1-benchmark" < run-search-hyperparameters.sh
# 12175707
# bsub -env CHARGE_METHOD="am1bcc",DATASET_NAME="openff-benchmark",FEATURE_NAME="v1-simple",HYPERPARAMETER_NAME="best-am1bcc-v1-benchmark" < run-search-hyperparameters.sh
# 12175708
# bsub -env CHARGE_METHOD="am1",DATASET_NAME="openff-benchmark",FEATURE_NAME="v5-espaloma-gasteiger",HYPERPARAMETER_NAME="best-am1-v5" < run-search-hyperparameters.sh
# 12184214
# bsub -env CHARGE_METHOD="am1bcc",DATASET_NAME="openff-benchmark",FEATURE_NAME="v5-espaloma-gasteiger",HYPERPARAMETER_NAME="best-am1bcc-v5" < run-search-hyperparameters.sh
# 12184215


# bsub -env CHARGE_METHOD="am1",DATASET_NAME="openff-benchmark",FEATURE_NAME="v6-espaloma-nores",HYPERPARAMETER_NAME="best-am1-v6-benchmark" < run-search-hyperparameters.sh
# 12202016
# bsub -env CHARGE_METHOD="am1bcc",DATASET_NAME="openff-benchmark",FEATURE_NAME="v6-espaloma-nores",HYPERPARAMETER_NAME="best-am1bcc-v6-benchmark" < run-search-hyperparameters.sh
# 12202019
# bsub -env CHARGE_METHOD="am1",DATASET_NAME="openff-benchmark",FEATURE_NAME="v7-espaloma-gasteiger-nores",HYPERPARAMETER_NAME="best-am1-v7-benchmark" < run-search-hyperparameters.sh
# 12202020
# bsub -env CHARGE_METHOD="am1bcc",DATASET_NAME="openff-benchmark",FEATURE_NAME="v7-espaloma-gasteiger-nores",HYPERPARAMETER_NAME="best-am1bcc-v7-benchmark" < run-search-hyperparameters.sh
# 12202021

