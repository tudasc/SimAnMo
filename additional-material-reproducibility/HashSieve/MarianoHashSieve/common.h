#ifndef COMMON_H
#define COMMON_H

//#######################

// N: Dimension
// SEED: Seed of the input challenge basis
// PROBE: Level of probing
//		PROBE = 0: Do not use probing (fastest but uses more space)
//		PROBE = 1: 1-level probing (factor ~2 slower, factor [1 + k/2] less space)
//		PROBE = 2: 2-level probing (factor ~4 slower, factor [1 + k/2 + k(k-1)/8] less space)
// K: "Hash length"
//		Theoretically: K = 0.22n
// T: Number of hash tables
// 		Theoretically (without probing): T = 2^(0.129n)
//		Theoretically (1-level probing): T = 2^(0.129n) / [1 + k/2]
//		Theoretically (2-level probing): T = 2^(0.129n) / [1 + k/2 + k(k-1)/8]
// dim 50: K = 11, T = 91
// dim 60: K = 13, T = 223, MVIS = 18000
// dim 70: K = 15, T = 549

/*#define N 40
#define K 9
#define TABLES 36*/

/*#define N 50
#define K 11
#define TABLES 87*/

/*#define N 60
#define K 13
#define TABLES 214*/

/*#define N 70
#define K 17
#define TABLES 523*/

/*#define N 80
#define K 18
#define TABLES 1278*/

#define BUCKETS (1 << (K - 1))

#endif
