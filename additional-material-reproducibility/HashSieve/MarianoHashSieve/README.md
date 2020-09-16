# HashSieve
This is an implementation of the algorithm presented in the paper: Parallel (probable) lock-free HashSieve: a practical sieving algorithm for the SVP.

## Folder Structure
* results: Contains Score-P measurements of several executions solving the SVP of lattices with many different dimensions.
* lattices: Contains lattice bases pre-reduced with BKZ, delta=0.99 and beta=34

## Build
The original implementation has parameters (Number of tables, bucket size, dimension) hardcoded in the source.
The script "build.sh" calculates and defines the parameters based on the dimension of the lattice HashSieve is working on.
This also means that, for HashSieve to solve the SVP on a lattice with a different dimension, HashSieve needs to be recompiled.
To build HashSieve:

	./build.sh {dim}

with dim being the dimension of the lattice, e.g.

	./build.sh 50

## Run
	./hashsieve lattice-file number-of-threads
	e.g. ./hashsieve lattices/svpchallengedim40seed0.txt 32
