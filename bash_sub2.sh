#!/bin/bash

set -e

#### proprecess data #######
# dos2unix
cd sub2/data/raw
dos2unix -q *

# trim data
cd ../trimmed_set
Rscript trim_data.r
Rscript create_subset1.0.r
Rscript create_top_feature_list.r

# cross-validation
cd ../cv_set
Rscript cv_partition.r

# normalization
cd ../normalization/raw
./preprocess_data.sh
##############################


#### run predictions ########

## 1.breast ################

## 1.1 rna (generic)
cd ../../../prediction/breast/rna
./run_prediction.sh 

## 1.2 individual
cd ../../../prediction/breast/individual
./run_prediction.sh # 1cv - 6parallel - 204m

# 1.2.1 individual_1000feature
cd ../../../prediction/breast/individual_1000feature
./run_prediction.sh # 1cv - 6parallel - 29m
# 1.2.2 individual_100feature
cd ../../../prediction/breast/individual_100feature
./run_prediction.sh # 1cv - 3parallel - 16m
# 1.2.3 individual_10feature
cd ../../../prediction/breast/individual_10feature
./run_prediction.sh # 1cv - 3parallel - 12m
# 1.2.4 individual_go_expression
cd ../../../prediction/breast/individual_go_expression
./run_prediction.sh # 1cv - 6parallel - 109m

## 1.3 individual-transplant
cd ../../../prediction/breast/individual_transplant
./run_prediction.sh # 1cv - 6parallel - 437m

# 1.3.1 individual_0.8sample
cd ../../../prediction/breast/individual_0.8sample
./run_prediction.sh # 1cv - 6parallel - 161m
# 1.3.2 individual_0.6sample
cd ../../../prediction/breast/individual_0.6sample
./run_prediction.sh # 1cv - 3parallel - 227m
# 1.3.3 individual_0.4sample
cd ../../../prediction/breast/individual_0.4sample
./run_prediction.sh # 1cv - 3parallel - 149m

## 1.4 stacking
cd ../final
Rscript stack_model.r

## 1.5 feature importance
cd ../feature_importance
./run_prediction.sh

#######################################


## 2.ova ################

## 2.1 rna (generic)
cd ../../../prediction/ova/rna
./run_prediction.sh

## 2.2 individual
cd ../../../prediction/ova/individual
./run_prediction.sh # 1cv - 6parallel - 129m

# 2.2.1 individual_1000feature
cd ../../../prediction/ova/individual_1000feature
./run_prediction.sh # 1cv - 3parallel - 51m
# 2.2.2 individual_100feature
cd ../../../prediction/ova/individual_100feature
./run_prediction.sh # 1cv - 3parallel - 45m
# 2.2.3 individual_10feature
cd ../../../prediction/ova/individual_10feature
./run_prediction.sh # 1cv - 3parallel - 36m
# 2.2.4 individual_go_expression
cd ../../../prediction/ova/individual_go_expression
./run_prediction.sh # 1cv - 3parallel - 187m

## 2.3 individual-transplant
cd ../../../prediction/ova/individual_transplant
./run_prediction.sh # 1cv - 6parallel - 425m

# 2.3.1 individual_0.8sample
cd ../../../prediction/ova/individual_0.8sample
./run_prediction.sh # 1cv - 6parallel - 101m
# 2.3.2 individual_0.6sample
cd ../../../prediction/ova/individual_0.6sample
./run_prediction.sh # 1cv - 3parallel - 143m
# 2.3.3 individual_0.4sample
cd ../../../prediction/ova/individual_0.4sample
./run_prediction.sh # 1cv - 3parallel - 98m

## 2.4 stacking
cd ../final
Rscript stack_model.r

## 2.5 feature importance
cd ../feature_importance
./run_prediction.sh
Rscript create_gene_list.r

######################################


