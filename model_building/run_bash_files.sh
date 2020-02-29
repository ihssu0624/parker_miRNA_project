#!/bin/bash

for i in {1..2}
do
	for j in {1..662}
		do
		   : 
		   # do whatever on $i
		   echo $i
		   echo $j  
		   echo "#!/bin/bash
		#SBATCH -p general
		#SBATCH -N 1
		#SBATCH --time=5:00:00
		#SBATCH --mem=50g
		#SBATCH -n 1
		#SBATCH --cpus-per-task=1
		#SBATCH --output=outfile_$i_$j.out

		module add r/3.6.0
		srun Rscript miRNA_model_build.R $i $j" >> run_R_models_${i}_${j}.sh
	done
done
