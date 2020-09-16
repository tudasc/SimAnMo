#ifndef SAMPLER_HEADER
#define SAMPLER_HEADER

#include <stdio.h>
#include <stdlib.h>

#include <math.h>

#include <omp.h>

#include <immintrin.h>
#include <xmmintrin.h>
#include <atomic>
#include "common.h"

typedef short MY_DATA;

struct NodeC {
	MY_DATA __attribute__((aligned(16))) data[N];
	//	unsigned short pos_in_bucket[1278];
	long norm;
	struct NodeC *next;
	char lock;
//  std::atomic_flag lock = ATOMIC_FLAG_INIT;
};

int SampleZ(double c, double s_square, double t_, unsigned int seed,
		drand48_data *private_status);

void SampleLF(struct NodeC *pLF, int n_, int m_, double t_, double **coef_,
		double *s_prime_square_, double **mu_, long **B_, int seed,
		drand48_data *private_status);

#endif
