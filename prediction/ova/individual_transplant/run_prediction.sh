#!/bin/bash

num=6

for j in {1..5} # 5-fold cv
do
	# prepare data
	Rscript prepare_transplant_data_in_parallel.r ${j} ${num}
	
	# run in parallel
	for i in `seq 1 ${num}`; # number of parallel
	do
		cd o${j}${i}/
		python fast_prediction.py &
		cd ..
	done
	wait
	
	# cloak & concatenate
	Rscript cloak_and_concatenate_transplant_pred.r ${j} ${num}
done

Rscript stack_model.r 

