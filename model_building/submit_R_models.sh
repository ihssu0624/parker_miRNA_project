#!/bin/bash

for i in {1..2}
do
	for j in {1..662}
		do
		  sbatch run_R_models_${i}_${j}.sh
	done
done
