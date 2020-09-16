#!/bin/sh
cd ~/../../../g6k
source ./activate
dim_array="94 98 104 108 114 118 124"
seed_array="0 684 559 629 192"
repeat=3
echo "numthreads: "$NUM_THREADS
export REDUCE=1
for dim in $dim_array;
do
	for seed in $seed_array;
	do
		if [[ $seed -eq 0 ]];then
			seed_=00
		else
			seed_=$seed
		fi

		if [[ $dim -lt 100 ]];then
			dim_="0$dim"
		else 
			dim_=$dim
		fi
		challenge="svpchallenge-dim-${dim_}-seed-${seed_}.txt"
		command="./svp_challenge.py --dim $dim --threads $NUM_THREADS --seed $seed --path /home/ex56yseb/Projects/crossing/NewSieving/Codes/g6k/svpchallenge/$challenge"
		for (( i=0; i < $repeat; ++i ))
		do		
			$command 
			#2>&1 | tee -a $tmp_output
		done
	done
done
