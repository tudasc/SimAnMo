#!/bin/sh
cd ~/../../../g6k
source ./activate
dim_array="120"
seed=0
repeat=20
export REDUCE=1
for dim in $dim_array;
do
	echo "numthreads: "$NUM_THREADS
	challenge="svpchallenge-dim-${dim}-seed-00.txt"
	command="./svp_challenge.py --dim $dim --threads $NUM_THREADS --seed $seed --path /home/ex56yseb/Projects/crossing/NewSieving/Codes/g6k/svpchallenge/$challenge"
	for (( i=0; i < $repeat; ++i ))
	do		
		$command 
		#2>&1 | tee -a $tmp_output
	done
done
