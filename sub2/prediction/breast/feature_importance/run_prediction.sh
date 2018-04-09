#!/bin/bash

num=12

# prepare data
Rscript prepare_data_in_parallel.r 0 ${num}

# run in parallel
for i in `seq 1 ${num}`;
do
	cd o${i}/
#	python fast_prediction.py &
	python rf_save.py &
	cd ..
done
wait

# collect all models
mkdir model
for i in `seq 1 ${num}`;
do
    mv o${i}/*model model/
done

# extract feature importance
mkdir importance
python extract_importance.py
perl concatenate_fi.pl list_breast_rna_subset1.0.txt

