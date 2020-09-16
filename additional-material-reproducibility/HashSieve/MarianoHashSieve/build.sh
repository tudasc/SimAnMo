#!/bin/bash
if [ $# -eq 0 ]
	then
    	echo "Supply dimension of lattice"
	exit
fi

module load ntl/11.3.2.bigmem

dir=$(dirname "$(readlink -f "$0")")

N=$1
T=$(python -c "print(int(round(2**(0.129*$N))))")
K=$(python -c "print(int(round(0.22*$N)))")
echo -e "Compiling HashSieve with: N=$N, K=$K, T=$T"
make DEFS="-DN=$N -DTABLES=$T -DK=$K"
