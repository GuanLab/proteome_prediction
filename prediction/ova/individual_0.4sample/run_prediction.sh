#!/bin/bash

# it take about 4 hours to run 1-fold cv in 6 parallel

num=3

for j in {1..5} # 5-fold cv
do
	# prepare data
	Rscript prepare_data_in_parallel.r ${j} ${num}
	
	# run in parallel
	for i in `seq 1 ${num}`; # number of parallel
	do
		cd o${j}${i}/
		python fast_prediction.py & 
		cd ..
	done
	wait
	
	# cloak & concatenate
	Rscript cloak_and_concatenate_pred.r ${j} ${num}
done
