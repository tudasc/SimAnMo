#include "sampler.h"
#include "common.h"

#include <ctype.h>
#include <math.h>
#include <time.h>

#include <getopt.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <fstream>
#include <iostream>

#include <NTL/LLL.h>
#include <NTL/ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/vec_double.h>

#include <immintrin.h>
#include <limits.h>
#include <smmintrin.h>
#include <xmmintrin.h>

#include <malloc.h>

#include <stdbool.h>

#include <omp.h>

#include <sstream>
#include <vector>

///#include <papi.h>

#ifdef USING_SCOREP
#include <scorep/SCOREP_User.h>
#endif

#define ERROR_RETURN(retval)                                                   \
  {                                                                            \
    fprintf(stderr, "Error %d %s:line %d: \n", retval, __FILE__, __LINE__);    \
    exit(retval);                                                              \
  }

#define NUMBER_OF_COUNTERS 4

//#include "bkz.h"

//#include <pthread.h>

//#define OPT1

#ifdef OPT1
#include <algorithm>
#endif

NTL_CLIENT

#define SSE
#define LOCKING
#define NTL

long iteration;
long collisions;
long no_vectors;

/* Each hash table bucket contains a list of pointers to vectors */
/* The length indicates the sentinel number of occupied entries */
/* The size indicates the sentinel size of the underlying array */
struct bucket {
	struct NodeC **pointers;
	unsigned int size;
	unsigned int length;
	char lock;
	struct bucket *next;
};

/* HashSieve-specific variables */
unsigned short AH[TABLES][K][2]; /* Hash matrices: only stores coordinates of non-zero stuff */
bucket HashTables[TABLES][BUCKETS]; /* Locality-sensitive hash tables */
unsigned long min_norm; /* Current minimum of squared vector norms */

// char HashTablesLOCKS[TABLES][BUCKETS];

unsigned seed_for_rand48 = 64;

// drand48_data status[32], *pts[32];

long maxI(long a, double b) {

	if (a > b)
		return a;

	return b;
}

/* Headers */
void add(struct NodeC *v, struct NodeC *w, long vw);
long ip(struct NodeC *v, struct NodeC *w);
void initHashes();
void MatIntFromMatZZ(long **A, const NTL::mat_ZZ B);
void MatZZFromMatInt(NTL::mat_ZZ B, long **A);
void MatDoubleFromMatRR(double **A, NTL::mat_RR B);
int max(int a, int b);
long lshash(struct NodeC *v, int t);
void bucket_add(bucket *b, struct NodeC *v);
void bucket_remove(bucket *b, struct NodeC *v /*, int vPos2*/);
void MatLongtoVect(long **A, struct NodeC *p, int no);
void VecZZToListPoint(const NTL::vec_ZZ &v, struct NodeC *p);
void convertFromStr(NTL::vec_ZZ &value, const std::string &s);
NTL::mat_ZZ randomizeMatrix(NTL::mat_ZZ A, NTL::vec_ZZ randoms, int seed,
		int DIM);
void insertSortedStack(struct NodeC **pvt_stack, struct NodeC *e,
		long *pvt_stack_size);
void HashSieve(int rows, int cols, double t_, long **B_,
		double *s_prime_square_, double **mu_, struct NodeC **basisvectors,
		struct NodeC *best);

unsigned long omp_get_thread_num_wrap() {

	unsigned long tid = omp_get_thread_num();

	return tid;
}

using namespace std;

int main(int argc, char **argv) {

	no_vectors = N;
	min_norm = LONG_MAX;

	MY_DATA data[N];

	NTL::mat_ZZ B;
	std::ifstream input_file(argv[1]);
	if (input_file.is_open()) {
		input_file >> B;
		input_file.close();
	} else {
		printf("File was not read\n");
		exit(-1);
	}

	double startBKZ = 0, endBKZ = 0;

	startBKZ = omp_get_wtime();

	int window;

	if (N > 74)
		window = 22;
	else
		window = 20;

	NTL::G_BKZ_FP(B, 0.99, window);

	endBKZ = omp_get_wtime();

	printf("Execution time of BKZ, with window = %d = %f.\n", window,
			endBKZ - startBKZ);

	initHashes();

	/* SAMPLER INITIALIZATION */

	NTL::mat_RR muq;
	NTL::vec_RR Bstar_squareq;
	NTL::ComputeGS(B, muq, Bstar_squareq);

	int rows = B.NumRows(), cols = B.NumCols(), i;

	long **B_ = (long **) malloc(B.NumCols() * sizeof(long *));

	for (int row = 0; row < rows; row++)
		B_[row] = (long *) malloc(B.NumRows() * sizeof(long));

	double **mu_ = (double **) malloc((B.NumCols()) * sizeof(double *));

	for (int row = 0; row < rows; ++row)
		mu_[row] = (double *) malloc(B.NumRows() * sizeof(double));

	// Datatype conversion:
	MatIntFromMatZZ(B_, B);

	MatDoubleFromMatRR(mu_, muq);

	double *s_prime_square_ = (double *) malloc(rows * (sizeof(double)));

	long max_star_sqr_norm = 0;

	for (i = 0; i < rows; ++i) {
		max_star_sqr_norm = max(max_star_sqr_norm, to_long(Bstar_squareq[i]));
	}

	double t_ = log(rows);

	if (N >= 80)
		t_ = log(rows) / 70;

	else
		t_ = log(rows) / 20;

	double s_square = max_star_sqr_norm * t_;

	for (i = 0; i < rows; ++i) {
		s_prime_square_[i] = s_square / to_double(Bstar_squareq[i]);
	}

	/* SAMPLER INITIALIZATION */

	iteration = 0;

	double start = 0, end = 0;

	int NTHREADS = atoi(argv[2]);

	omp_set_num_threads(NTHREADS);

	if (argc < 5)
		printf(
				"\t Lattice %s in dimention %d. K = %d, TABLES = %d and Buckets = %d:\n\n",
				argv[1], N, K, TABLES, BUCKETS);
	else
		printf("\t Lattice %s in dimention %d, with target norm %d:\n\n",
				argv[1], N, atoi(argv[4]));

	// Pointer to the best vector at any instant
	struct NodeC *shortest = (NodeC *) malloc(sizeof(struct NodeC));
	shortest->next = NULL;
	shortest->lock = 0;

	struct NodeC **basisvectors = (struct NodeC **) malloc(
			N * sizeof(struct NodeC *));
	for (int v = 0; v < N; v++) {
		basisvectors[v] = (struct NodeC *) malloc(sizeof(struct NodeC));
		basisvectors[v]->next = NULL;
		basisvectors[v]->lock = 0;
	}

	/* PAPI Benchmarking */
	/*
	 int retval, N_THREADS = 32;

	 if ((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT)
	 ERROR_RETURN(retval);
	 */
	//if (
	/*
	 PAPI_thread_init(
	 (unsigned long (*)(void))(omp_get_num_threads))         !=
	 PAPI_OK ) printf("ERROR WITH PAPI_THREADS\n");
	 */
	//  PAPI_thread_init(omp_get_thread_num_wrap) != PAPI_OK)
	// printf("ERROR WITH PAPI_THREADS\n");
	// printf("PAPI sparked off with %ld threads.\n",t_p);
	/*
	 long long *values =
	 (long long *)malloc(sizeof(long long) * N_THREADS * NUMBER_OF_COUNTERS);
	 int *EventSet = (int *)malloc(sizeof(int) * N_THREADS);
	 memset(values, 0, sizeof(long long) * N_THREADS * NUMBER_OF_COUNTERS);

	 #pragma omp parallel private(retval)
	 {
	 int tid = omp_get_thread_num();
	 EventSet[tid] = PAPI_NULL;
	 if ((retval = PAPI_create_eventset(&EventSet[tid])) != PAPI_OK)
	 ERROR_RETURN(retval);
	 if ((retval = PAPI_add_event(EventSet[tid], PAPI_L2_TCM)) != PAPI_OK)
	 ERROR_RETURN(retval);
	 if ((retval = PAPI_add_event(EventSet[tid], PAPI_L2_TCA)) != PAPI_OK)
	 ERROR_RETURN(retval);
	 if ((retval = PAPI_add_event(EventSet[tid], PAPI_L3_TCM)) != PAPI_OK)
	 ERROR_RETURN(retval);
	 if ((retval = PAPI_add_event(EventSet[tid], PAPI_L3_TCA)) != PAPI_OK)
	 ERROR_RETURN(retval);
	 if ((retval = PAPI_start(EventSet[tid])) != PAPI_OK)
	 ERROR_RETURN(retval);
	 }*/

	/* HashSieve call */

	start = omp_get_wtime();

	HashSieve(rows, cols, t_, B_, s_prime_square_, mu_, basisvectors, shortest);

	end = omp_get_wtime();

// stop papi counters and get values
	/*#pragma omp parallel private(retval)
	 {
	 int tid = omp_get_thread_num();
	 if ((retval = PAPI_stop(EventSet[tid],
	 &values[tid * NUMBER_OF_COUNTERS])) != PAPI_OK)
	 ERROR_RETURN(retval);
	 if ((retval = PAPI_reset(EventSet[tid])) != PAPI_OK)
	 ERROR_RETURN(retval);
	 }

	 // this stores the values of each counter for each thread, in the values
	 // array. to obtain aggregate results we need to sum up everything
	 long long sum[NUMBER_OF_COUNTERS];

	 for (i = 0; i < NUMBER_OF_COUNTERS; i++)
	 sum[i] = 0;

	 for (i = 0; i < N_THREADS; i++) {
	 sum[0] += values[i * NUMBER_OF_COUNTERS + 0];
	 sum[1] += values[i * NUMBER_OF_COUNTERS + 1];
	 sum[2] += values[i * NUMBER_OF_COUNTERS + 2];
	 sum[3] += values[i * NUMBER_OF_COUNTERS + 3];
	 }

	 printf(" PAPI TIMINGS: \n\n\t PAPI_L2_TCM = %ld\n\t PAPI_L2_TCA = %ld\n\t "
	 "PAPI_L3_TCM = %ld\n\t PAPI_L3_TCA = %ld \n\n\n",
	 sum[0], sum[1], sum[2], sum[3]);

	 // cleanup
	 #pragma omp parallel private(retval)
	 {
	 int tid = omp_get_thread_num();
	 if ((retval = PAPI_remove_event(EventSet[tid], PAPI_L2_TCM)) != PAPI_OK)
	 ERROR_RETURN(retval);
	 if ((retval = PAPI_remove_event(EventSet[tid], PAPI_L2_TCA)) != PAPI_OK)
	 ERROR_RETURN(retval);
	 if ((retval = PAPI_remove_event(EventSet[tid], PAPI_L3_TCM)) != PAPI_OK)
	 ERROR_RETURN(retval);
	 if ((retval = PAPI_remove_event(EventSet[tid], PAPI_L3_TCA)) != PAPI_OK)
	 ERROR_RETURN(retval);
	 if ((retval = PAPI_destroy_eventset(&(EventSet[tid]))) != PAPI_OK)
	 ERROR_RETURN(retval);
	 }
	 */
	/* Prints solution */

	printf("\t ******** SOLUTION ****** Vector norm  = %ld!\n\n[", min_norm);

	for (i = 0; i < N; i++)
		if (i == 0 || i == N - 1)
			printf("%d", shortest->data[i]);
		else
			printf(" %d ", shortest->data[i]);

	printf("]\n\n\tAlgorithm ended in %ld iterations and with %ld collisions. "
			"No. vectors = %ld.\n", iteration, collisions, no_vectors);
	printf("\n\t\tHashSievet took (N=%d):%f seconds.\n\n", N, end - start);

	/* Provides process info */

	char str[80];

	// sprintd(str, "pmap -x %d\0", getpid());
	// system(str);

	sprintf(str, "ps -v %d\0", getpid());
	system(str);

	/*
	 FILE *f = fopen("density.txt", "w");



	 // Density map Hash tables
	 for(int i=0;i<TABLES;i++){
	 for(int j=0;j<BUCKETS;j++){
	 fprintf(f,"%d ",HashTables[i][j].length);
	 }
	 fprintf(f,"\n");
	 }


	 fclose(f);
	 */

	exit(1);

	// return 0;
}

void HashSieve(int rows, int cols, double t_, long **B_,
		double *s_prime_square_, double **mu_, struct NodeC **basisvectors,
		struct NodeC *shortest) {

#ifdef USING_SCOREP
	SCOREP_USER_REGION( "HashSieve", SCOREP_USER_REGION_TYPE_FUNCTION )
#endif
	// Auxilitiary pointer to the best vector in the algorithm at every instant
	struct NodeC *best;

	// Marks the end of the execution (the thread that reaches the no. of
	// collisions first notifies the others with this variable)
	int done = 0;

#pragma omp parallel shared(done, rows, cols, s_prime_square_, mu_, B_)
	{

		// Structure necessary for the sampler
		double *coef_ = (double *) malloc(rows * (sizeof(double)));

		/* Declares pool of vectors to use for samples */

		struct NodeC *pool = (struct NodeC *) malloc(
				100000 * sizeof(struct NodeC));

		for (int u = 0; u < 100000; u++) {
			pool[u].lock = 0;
			pool[u].next = NULL;
		}

		// "Pointer" to the position of the pool
		int pos_pool = 0;

		// Private pointer to the current sample (simplifies the process)
		struct NodeC *v;

		// Thread ID
		int tid = omp_get_thread_num();

		/* BEGIN: Initialization of private states for drand48_r() */

		drand48_data status, *pts;

		pts = &status;
		srand48_r(seed_for_rand48 * tid, pts);

		int private_iterations = 0;

		double start_rand, end_rand;

		int retLR = 1;

		int maxSizeAux = 0;

		/* END: Initialization of private states for rand48() */

		// Pointer to the begining of the (private) stack
		struct NodeC *pvt_stack = NULL;

		// Current size of the stack
		long pvt_stack_size = 0;

		struct NodeC aux_one_original;
		struct NodeC aux_one_original2;
		struct NodeC aux_two_minimal;

		/* Load basis vectors to the private stack of each thread */

#pragma omp for schedule(static) nowait
		for (int i = 0; i < N; ++i) {
			// Convert one basis vector to one "real" vector
			MatLongtoVect(B_, basisvectors[i], i);

			// Move basisvectors[i] onto the stack:
			// Inserted ordered:
			//      insertSortedStack(&pvt_stack, basisvectors[i],
			//&pvt_stack_size);
			// normal insertion
			/*
			 cout << "[";
			 for(int j=0;j<N;j++){
			 if(j!=N-1) cout <<
			 basisvectors[i]->data[j] << " " ;
			 else cout << basisvectors[i]->data[j] <<
			 "]" << endl ;
			 }
			 */

			basisvectors[i]->next = pvt_stack;
			pvt_stack = basisvectors[i];
			pvt_stack_size++;
		}

		/* Begin: HashSieve */

		long priv_cols = 0;
		long priv_num_vectors = 0;
		// Auxiliary norm for repeated samples of zero vectors
		long aux_norm_null = 0;

		int no_candidates;
		struct NodeC **candidates;
		struct NodeC *w;
		long vw;
		unsigned long vwAbs;
		unsigned long vwAbsDouble;

		int t, tt, valHash, j, it, sampled = 0, batch_cols = 0;

		long index;

		/* Main iteration loop */

#pragma omp flush(done)
		while (!done) {

			newit:
			// No. of iterations done by one thread - relevant for statistic purposes
			// only
			private_iterations++;

			/* SAMPLE A VECTOR OR PICK ONE VECTOR FROM THE PRIVATE STACK */

			// If there are no vectors in the stack:
			// if(stack_pool_pos==0){
			if (pvt_stack_size == 0) {
				// if(stackSizeLF==0){

				/*
				 if(sampled>0){
				 v = &pool[ pos_pool-sampled ];
				 sampled--;
				 priv_num_vectors++;
				 }
				 else{

				 #pragma omp atomic
				 collisions+=batch_cols;

				 batch_cols=0;

				 //cout << "A new batch will be executed" << endl;
				 for(it=0;it<500;it++){
				 v = &pool[ pos_pool ];
				 pos_pool++;
				 sampled++;
				 SampleLF(v,rows,cols,t_,&coef_,s_prime_square_,mu_,B_,tid,pts);
				 aux_norm_null = v->norm;
				 while(aux_norm_null==0){
				 //#pragma omp atomic
				 //collisions++;
				 batch_cols++;
				 SampleLF(v,rows,cols,t_,&coef_,s_prime_square_,mu_,B_,tid,pts);
				 aux_norm_null = v->norm;
				 }
				 }
				 //cout << "The batch was executed" << endl;
				 v = &pool[ pos_pool-sampled ];
				 sampled--;
				 priv_num_vectors++;

				 }
				 */

				// We use another vector from our pool (a malloc here is too expensive)
				v = &pool[pos_pool];
				pos_pool++;

				// No. of vectors in the entire system (priv_num_vectors is the no. of
				// vectors in one thread)
				priv_num_vectors++;

				// Sample one vector (v is the handler, the vector is stored in pool[
				// pos_pool ]
				SampleLF(v, rows, cols, t_, &coef_, s_prime_square_, mu_, B_,
						tid, pts);

				// Re-sample the vector until the norm is bigger than 0:
				aux_norm_null = v->norm;
				while (aux_norm_null == 0) {
#pragma omp atomic
					collisions++;
					SampleLF(v, rows, cols, t_, &coef_, s_prime_square_, mu_,
							B_, tid, pts);
					aux_norm_null = v->norm;
				}

				// ONLY FOR DEBUG
				/*
				 for(int i=0; i<N ;i++)
				 cout << v->data[i] << " ";

				 index = lshash(v, 5);

				 cout << "index = " << index << endl;
				 getchar();
				 */
				//
			} else {
				// Grab one vector from the stack:
				v = pvt_stack;
				pvt_stack = pvt_stack->next;
				pvt_stack_size--;
			}

			// ========================================================

			/* Here we check each table for candidate near list vectors
			 That means that we calculate the hash of our sample (either
			 from scratch or popped from the stack and we go through all
			 the hashtables to reduce it and use it to reduce the bucket
			 vectors, i.e. the vectors in the buckets of the hashtables. */

			// ========================================================
			rep:
			// ========================================================
			/* We copy v to two different auxiliary vectors (aux_one_
			 original and aux_one_original2) that will be used to record
			 the minimum representative of v when traversing the entire
			 bucket even after reducing v once inside one bucket. This is
			 noted below as LOCALITY_PROCEDURE;           */
			// ========================================================
			memcpy(aux_one_original.data, v->data, N * sizeof(MY_DATA));
			memcpy(aux_one_original2.data, v->data, N * sizeof(MY_DATA));
			aux_one_original.norm = v->norm;
			aux_one_original2.norm = v->norm;
			aux_one_original.lock = 0;
			// aux_two_minimal.norm = LONG_MAX;
			aux_two_minimal.lock = 0;

			// For each hash table:
			for (t = 0; t < TABLES; t++) {

				// Calculate the hash value of in hash table t:
				index = lshash(v, t);

				// For each candidate vector:
				for (j = 0; j < HashTables[t][index].length; j++) {

					/*
					 * If this bucket is being used (element added or removed),
					 * we block till we get its lock;
					 */

#ifdef LOCKING
					while (__builtin_expect(
							!(__sync_bool_compare_and_swap(
									&(HashTables[t][index].lock), 0, 1)), 0))
						;
#endif
					w = HashTables[t][index].pointers[j];
#ifdef LOCKING
					HashTables[t][index].lock = 0;
#endif

					/*
					 * Try to get lock of w. We skip w if no lock is gotten
					 * but we do not release v, as we want to have it locked
					 * in the next iteration of j; if j moves forward with
					 * continue in every iteration, it doesn't matter because
					 * v is released after this look anyway.
					 */
#ifdef LOCKING
					if (__builtin_expect(
							!(__sync_bool_compare_and_swap(&(w->lock), 0, 1)),
							0)) {
						continue;
					}
#endif

					// Calculate the inner product between v and w:
					vw = ip(v, w);

					// FOR PROFILER REASONS (if we activate this the code runs slow like a
					// dog)
					//#pragma omp atomic
					// innerps++;

					vwAbs = (vw > 0 ? vw : -vw);
					vwAbsDouble = (vwAbs << 1);

					if (vwAbsDouble > w->norm) {

						add(v, w, vw);

						/*
						 * If v is modified, we release w's lock, and one jump
						 * to iterate over the next hashtable (keeping v's lock);
						 */
#ifdef LOCKING
						w->lock = 0;
#endif
						if (v->norm == 0) {
#pragma omp atomic
							collisions++;

							priv_num_vectors--;

							goto newit;
						}

						/* =================================================================
						 */
						/* ======================= LOCALITY_PROCEDURE ======================
						 */
						/* =================================================================
						 */

						/*
						 * The original algorithm would jump at this point, because v has
						 * been
						 * changed. However, in this version we look at the whole bucket,
						 * always
						 * keeping the shortest version of v; This is because (1) we might
						 * encounter better versions of v, (2) we could reduce more
						 * candidate
						 * vectors (w) against v and (3) this favors locality because the
						 * bucket
						 * is loaded to cache.
						 *
						 * Whenever we find a better version of v, we copy it to
						 * aux_two_mininal,
						 * wherein we store its shortest version;
						 */

						/*
						 * NOTE: The whole locality procedure can be commented out and no
						 * further
						 * changes are necessary for the code to work;
						 */

						// Calculate the new hash value of v, which changed because v
						// changed:
						long tmpvhash = lshash(v, 0);

						// Prefetch that memory position for the first hash table:
						__builtin_prefetch(&(HashTables[0][tmpvhash]), 1, 0);

						// Save in "aux_two_minimal" the current minimum of v:
						memcpy(aux_two_minimal.data, v->data,
								N * sizeof(MY_DATA));
						aux_two_minimal.norm = v->norm;

						// For the remaining candidate vectors:
						for (int z = j + 1; z < HashTables[t][index].length;
								z++) {
#ifdef LOCKING
							while (!(__sync_bool_compare_and_swap(
									&(HashTables[t][index].lock), 0, 1)))
								;
#endif
							w = HashTables[t][index].pointers[z];
#ifdef LOCKING
							HashTables[t][index].lock = 0;
#endif

#ifdef LOCKING
							if (!(__sync_bool_compare_and_swap(&(w->lock), 0, 1))) {
								continue;
							}
#endif

							// Calculate the inner product between the original v and the
							// candidate vector w:
							vw = ip(&aux_one_original, w);

							vwAbs = (vw > 0 ? vw : -vw);
							vwAbsDouble = (vwAbs << 1);

							// If vwAbsDouble > w->norm holds than v will get smaller
							if (vwAbsDouble > w->norm) {

								/*
								 * We could calculate the resultant norm of aux_one_original
								 * without actually
								 * modifying the vector ifself (i.e. call function "add") with:
								 * if(vw > 0) norm_n = aux_one_original.norm + w->norm - 2 * vw;
								 * else norm_n = aux_one_original.norm + w->norm + 2 * vw;
								 * but this does not deliver performance gains, and so the most
								 * readable version
								 * was kept instead;
								 */

								add(&aux_one_original, w, vw);

								// If v will get smaller, than we actually reduce the original v
								// (in aux_one_original) and store it in aux_two_minimal
								if (aux_one_original.norm
										< aux_two_minimal.norm) {
									memcpy(aux_two_minimal.data,
											aux_one_original.data,
											N * sizeof(MY_DATA));
									aux_two_minimal.norm =
											aux_one_original.norm;
								}

								// We also copy the original v to "aux_one_original", which is
								// used for actual reductions:
								memcpy(aux_one_original.data,
										aux_one_original2.data,
										N * sizeof(MY_DATA));
								aux_one_original.norm = aux_one_original2.norm;
							}

							// The bucket vector w is to be reduced against (the original) v:
							if (vwAbsDouble > aux_one_original.norm) {

								long hashesw[TABLES], wHash;

								hashesw[0] = lshash(w, 0);

								__builtin_prefetch(&(HashTables[0][hashesw[0]]),
										1, 0);

								for (int uu = 1; uu < TABLES; uu++)
									hashesw[uu] = lshash(w, uu);

								for (tt = 0; tt < TABLES - 1; tt++) {
									wHash = hashesw[tt];

									__builtin_prefetch(
											&(HashTables[tt + 1][hashesw[tt + 1]]),
											1, 0);

									while (!((__sync_bool_compare_and_swap(
											&(HashTables[tt][wHash].lock), 0, 1)))) {
									}

									bucket_remove(&HashTables[tt][wHash], w);

									HashTables[tt][wHash].lock = 0;
								}

								wHash = hashesw[tt];

								while (!((__sync_bool_compare_and_swap(
										&(HashTables[tt][wHash].lock), 0, 1)))) {
								}

								bucket_remove(&HashTables[tt][wHash], w);

								HashTables[tt][wHash].lock = 0;

								// Reduce w against (the original) v
								add(w, &aux_one_original, vw);

								if (w->norm > 0) {
									w->next = pvt_stack;
									pvt_stack = w;
									pvt_stack_size++;
								} else {
#pragma omp atomic
									collisions++;

									priv_num_vectors--;
								}
							}

							w->lock = 0;
						}

						// Save the minimum represented version of v in v:
						memcpy(v->data, aux_two_minimal.data,
								N * sizeof(MY_DATA));
						v->norm = aux_two_minimal.norm;

						/* =================================================================
						 */
						/* ================= END OF LOCALITY_PROCEDURE =====================
						 */
						/* =================================================================
						 */

						goto rep;
					}

					/* If the code didn't enter the previous if statment,
					 then we will eventually use the next candidate. Thus,
					 we can prefetch its pointer and data. In principle,
					 its pointer is fetched by the compiler automatically,
					 but its data of the vector ain't. However, activating
					 the prefetch below does not increase performance. */

					//__builtin_prefetch
					//(&(HashTables[t][index].pointers[j+1]->data),1,1);
					// Candidate w is to be reduced:
					if (vwAbsDouble > v->norm) {

						int hashesw[TABLES];

						if (t == 0)
							hashesw[0] = index;
						else
							hashesw[0] = lshash(w, 0);

						__builtin_prefetch(&(HashTables[0][hashesw[0]]), 1, 0);

						// Calculate the hashvalue of w in tt
						for (int uu = 1; uu < TABLES; uu++)
							if (uu == t)
								hashesw[uu] = index;
							else
								hashesw[uu] = lshash(w, uu);

						// This for goes originally until TABLES, but it is here modified (goes
						// till TABLES-1) due to prefetching
						for (tt = 0; tt < TABLES - 1; tt++) {

							valHash = hashesw[tt];

							// Prefetch the next hash table in the loop
							__builtin_prefetch(
									&(HashTables[tt + 1][hashesw[tt + 1]]), 1,
									0);
#ifdef LOCKING
							while (!((__sync_bool_compare_and_swap(
									&(HashTables[tt][valHash].lock), 0, 1)))) { /*cout << "waiting...2" << endl;*/
							}
#endif
							bucket_remove(&HashTables[tt][valHash], w);
#ifdef LOCKING
							HashTables[tt][valHash].lock = 0;
#endif
						}

						valHash = hashesw[tt];

// Last iteration of the previous loop:
#ifdef LOCKING
						while (!((__sync_bool_compare_and_swap(
								&(HashTables[tt][valHash].lock), 0, 1)))) { /*cout << "waiting...2" << endl;*/
						}
#endif
						bucket_remove(&HashTables[tt][valHash], w);
#ifdef LOCKING
						HashTables[tt][valHash].lock = 0;
#endif

						add(w, v, vw);

						if (w->norm > 0) {
							//              insertSortedStack(&pvt_stack,
							//w, &pvt_stack_size);
							w->next = pvt_stack;
							pvt_stack = w;
							pvt_stack_size++;
						} else {
#pragma omp atomic
							collisions++;

							priv_num_vectors--;
						}
					}

					/*
					 * The candidate list overs here, and the lock of w
					 * has to be released;
					 */
#ifdef LOCKING
					w->lock = 0;
#endif

					/* End of candidate loop (for) */
				}

// One layer of probing: visiting all adjacent buckets as well
#ifdef PROBE
				for (int bit = 0; bit < K; bit++) {

					// Compute v's hash value

					int vHashp = lshash(v, t) ^ (1 << bit);

					if (vHashp >= BUCKETS)
					vHashp = 2 * BUCKETS - vHashp - 1;

					for (j = 0; j < HashTables[t][vHashp].length; j++) {

						while (!(__sync_bool_compare_and_swap(&(HashTables[t][vHashp].lock),
												0, 1)))
						;

						w = HashTables[t][vHashp].pointers[j];

						HashTables[t][vHashp].lock = 0;
#ifdef LOCKING
						if (!(__sync_bool_compare_and_swap(&(w->lock), 0, 1))) {
							continue;
						}
#endif
						// Go through the list to find reducing vectors

						vw = ip(v, w);

						//#pragma omp atomic
						// innerps++;

						vwAbs = (vw > 0 ? vw : -vw);

						vwAbsDouble = (vwAbs << 1);

						// Reduce v with w if possible
						if (vwAbsDouble > w->norm) {

							add(v, w, vw);
#ifdef LOCKING
							w->lock = 0;
#endif
							if (v->norm == 0) {
#pragma omp atomic
								collisions++;

								priv_num_vectors--;

								goto newit;
							}

							goto rep;
						}

						int wHash;

						// Reduce w with v if possible
						if (vwAbsDouble > v->norm) {

							// Remove w from the hash tables
							// for(int tt = 0; tt < TABLES; tt++){

							long hashesw[TABLES];

							hashesw[0] = lshash(w, 0);

							//__builtin_prefetch (&(HashTablesLOCKS[0][hashesw[0]]),1,0);
							__builtin_prefetch(&(HashTables[0][hashesw[0]]), 1, 0);

							// Calculate the hashvalue of w in tt (this might have been
							// calculated above, but its better to re-calc)
							// because accessing more memory at this point will be painful
							for (int uu = 1; uu < TABLES; uu++)
							hashesw[uu] = lshash(w, uu);

							for (tt = 0; tt < TABLES - 1; tt++) {
								wHash = hashesw[tt];

								//__builtin_prefetch
								//(&(HashTablesLOCKS[tt+1][hashesw[tt+1]]),1,0);
								__builtin_prefetch(&(HashTables[tt + 1][hashesw[tt + 1]]), 1,
										0);

								while (!((__sync_bool_compare_and_swap(
																&(HashTables[tt][wHash].lock), 0, 1)))) {
								}
								bucket_remove(&HashTables[tt][wHash],
										w); //, w->pos_in_bucket[tt]);
								HashTables[tt][wHash].lock = 0;
							}

							wHash = hashesw[tt];

							while (!((__sync_bool_compare_and_swap(
															&(HashTables[tt][wHash].lock), 0, 1)))) {
							}
							bucket_remove(&HashTables[tt][wHash],
									w); //, w->pos_in_bucket[tt]);
							HashTables[tt][wHash].lock = 0;

							//  wHash = lshash(w, tt);

							//  while(!((__sync_bool_compare_and_swap(&(HashTables[tt][wHash].lock),0,1)))){
							///*cout << "waiting... 2" << endl;*/ }

							//  bucket_remove(&HashTables[tt][wHash], w);

							//  HashTables[tt][wHash].lock = 0;

							//}

							// Reduce w with v
							add(w, v, vw);

							if (w->norm > 0) {
								w->next = pvt_stack;
								pvt_stack = w;
								pvt_stack_size++;
								//                insertSortedStack(&pvt_stack,
								//w, &pvt_stack_size);
							} else {
#pragma omp atomic
								collisions++;

								priv_num_vectors--;
							}
						}
#ifdef LOCKING
						w->lock = 0;
#endif
					}
				}
#endif
				/* END: PROBING */
			} // End loop hashtables

			/*
			 * We lock v here since it will be visible by other threads
			 * at this point, because it will be added to the Hashtables
			 */

#ifdef LOCKING
			while (!(__sync_bool_compare_and_swap(&(v->lock), 0, 1))) {
				cout << "locking" << endl;
			}
#endif

			// __builtin_prefetch (&no_vectors,1,3);

			/* Since all the hash values of v will be calculated in the
			 next "t" loop above, we can calculate all the hashes upfront
			 and prefetch the Hash table positions and locks so they are
			 in cache when the right iteration comes; This improves the
			 code in about 10%. */

			int hashesv[TABLES];

			hashesv[0] = lshash(v, 0);

			//__builtin_prefetch (&(HashTablesLOCKS[0][hashesv[0]]),1,0);
			__builtin_prefetch(&(HashTables[0][hashesv[0]]), 0, 0);

			for (int uu = 1; uu < TABLES; uu++) {

				// if(uu!=TABLES-1){
				//  __builtin_prefetch (&(AH[uu+1][0][0]),0,0);
				//  __builtin_prefetch (&(AH[uu+1][0][1]),0,0);
				//}

				hashesv[uu] = lshash(v, uu);
			}

			// This loop is originaly executed from 0 to TABLES but runs till TABLES-1 due to
			// prefetching
			for (t = 0; t < TABLES - 1; t++) {

				valHash = hashesv[t];

				//__builtin_prefetch (&(HashTablesLOCKS[t+1][hashesv[t+1]]),1,0);
				__builtin_prefetch(&(HashTables[t + 1][hashesv[t + 1]]), 0, 0);

#ifdef LOCKING
				// One has to get the lock of the hashtable and bucket were v is to be
				// inserted
				// while(!((__sync_bool_compare_and_swap(&(HashTablesLOCKS[t][valHash]),0,1)))){
				// /* cout << "waiting... 3" << endl;*/  }
				while (!((__sync_bool_compare_and_swap(
						&(HashTables[t][valHash].lock), 0, 1)))) { /* cout << "waiting... 3" << endl;*/
				}
#endif
				// Add v to the right bucket in Hash table t
				/*v->pos_in_bucket[t] =*/bucket_add(&HashTables[t][valHash], v);
#ifdef LOCKING
				// Release the lock
				// HashTablesLOCKS[t][valHash] = 0;
				HashTables[t][valHash].lock = 0;
#endif
			}

			valHash = hashesv[t];

#ifdef LOCKING
			// One has to get the lock of the hashtable and bucket were v is to be
			// inserted
			// while(!((__sync_bool_compare_and_swap(&(HashTablesLOCKS[t][valHash]),0,1)))){
			// /* cout << "waiting... 3" << endl;*/  }
			while (!((__sync_bool_compare_and_swap(
					&(HashTables[t][valHash].lock), 0, 1)))) { /* cout << "waiting... 3" << endl;*/
			}
#endif
			// Add v to the right bucket in Hash table t
			/*v->pos_in_bucket[t] =*/bucket_add(&HashTables[t][valHash], v);
#ifdef LOCKING
			// Release the lock
			// HashTablesLOCKS[t][valHash] = 0;
			HashTables[t][valHash].lock = 0;
#endif

			if (v->norm < min_norm) {
				best = v;
				min_norm = v->norm;
				cout << "New best vector found, with norm = " << min_norm
						<< endl;
				std::cout << "\t no_vectors = " << no_vectors << std::endl;
				std::cout << "\t Collisions = " << collisions << std::endl;
			}

			/*
			 * Release lock of v, since v is no longer used
			 */
#ifdef LOCKING
			v->lock = 0;
#endif

#pragma omp atomic
			no_vectors += priv_num_vectors;

			priv_num_vectors = 0;

			if (collisions
					>= 0.1 * no_vectors + 200 /*|| min_norm <= 5162838*/) {

				done = 1;
#pragma omp flush(done)
			}

			/* End of while(done) */
		}

#pragma omp atomic
		iteration += private_iterations;
		private_iterations = 0;
	}

	memcpy(shortest->data, best->data, N * sizeof(MY_DATA));
}

// Initialize the hash vectors and the hash tables
void initHashes() {

#ifdef OPT1
	//
	int range[N];
	for (int i = 0; i < N; i++) {
		range[i] = i;
	}
#endif

	// Initialize hash tables as empty
	for (int t = 0; t < TABLES; t++) {

#ifdef OPT1
		//
		std::random_shuffle(range, range + N);
		for (int k = 0; k < K; k++) {
			AH[t][k][0] = range[2 * k];
			AH[t][k][1] = range[2 * k + 1];
		}
#else
		// Initialize random sparse hash vectors by choosing two non-zero entries
		for (int k = 0; k < K; k++) {
			AH[t][k][0] = (rand() % N);
			AH[t][k][1] = (rand() % N);
			while (AH[t][k][1] == AH[t][k][0])
				AH[t][k][1] = (rand() % N);
		}
#endif
		// Initialize empty hash buckets
		for (int b = 0; b < BUCKETS; b++) {
			HashTables[t][b].length = 0;
#ifndef LIST
			HashTables[t][b].size = 20;
#endif

#ifdef LIST
			HashTables[t][b].chaining = NULL;
#else
			HashTables[t][b].pointers = (struct NodeC **) malloc(
					20 * sizeof(struct NodeC *));
#endif
			//      HashTables[t][b].chaining = NULL;
			// HashTablesLOCKS[t][b] = 0;
			HashTables[t][b].lock = 0;
			HashTables[t][b].next = NULL;
		}
	}
}

void add(struct NodeC *v, struct NodeC *w, long vw) {

#ifdef SSE
	if (vw > 0) {

		int q = (vw + w->norm / 2) / w->norm;

		int16_t * const pA = v->data;
		int16_t * const pB = w->data;

#pragma ivdep
		for (int i = 0; i < N; i++)
			pA[i] -= q * pB[i];

		v->norm = v->norm + q * q * w->norm - 2 * q * vw;
	} else {
		register int i = 0;

		const int loopBound = N - 7;

		__m128i vsum, vecPi, vecCi, vecQCi;

		int16_t * const pA = v->data;
		int16_t * const pB = w->data;

		for (; i < loopBound; i += 8) {
			vecPi = _mm_load_si128((__m128i *) &(pA)[i]);
			vecCi = _mm_load_si128((__m128i *) &(pB)[i]);
			vecQCi = _mm_add_epi16(vecPi, vecCi);
			_mm_store_si128((__m128i *) &(pA)[i], vecQCi);
		}

		for (; i < N; i++)
			v->data[i] += pB[i];

		v->norm += w->norm + 2 * vw;
	}

#else
	if (vw > 0) {
		for (int i = 0; i < N; i++)
		v->data[i] -= w->data[i];
		v->norm += w->norm - 2 * vw;
	} else {
		for (int i = 0; i < N; i++)
		v->data[i] += w->data[i];
		v->norm += w->norm + 2 * vw;
	}
#endif
}

void bucket_add(bucket *b, struct NodeC *v) {

	if (b->length == b->size) {

		b->size = b->size + 20;
		struct NodeC **NewPointers;
		struct NodeC **tmp;

		NewPointers = (struct NodeC **) realloc(b->pointers,
				b->size * sizeof(struct NodeC *));

		b->pointers = NewPointers;
	}

	/*
	 * SORTED INSERSION; The bucket is traversed till a vector whose norm
	 * is bigger than than v's; The remaining vectors are shifted by one
	 * position and v is inserted in between;
	 */
	/*
	 if(b->length==0){

	 b->pointers[0] = v;
	 b->length++;

	 }
	 else{

	 int i = 0;

	 // This will traverse the bucket till a vector w s.t. w->norm
	 > v->norm is found
	 while(i < b->length){
	 if(v->norm > b->pointers[i]->norm)
	 i++;
	 else
	 break;
	 }

	 // Shifts the rest of the array by one position
	 int j = b->length;

	 while( j > i ){

	 b->pointers[j] = b->pointers[j-1];
	 j--;

	 }

	 // Inserts v in between
	 b->pointers[i] = v;
	 b->length++;

	 }
	 */

	/*
	 * UNSORTED INSERTION:
	 *
	 */

	b->pointers[b->length] = v;
	b->length++;
}

void bucket_remove(bucket *b, struct NodeC *v /*, int vPos2*/) {

	/*
	 if( b->pointers[vPos2] == v && vPos2 < b->length ){
	 b->length--;
	 b->pointers[vPos2] = b->pointers[b->length];
	 //    return vPos2;
	 }
	 else{
	 */
	int vPos = 0;
	/*
	 struct NodeC* sentinel = b->pointers[0];

	 if(sentinel->norm == v->norm){
	 b->pointers[0] = sentinel->next;
	 return;
	 }
	 if(!sentinel->next){ perror("Vector not found in bucket...\n");
	 return; }

	 while(sentinel->next->norm!=v->norm && sentinel->next)
	 sentinel = sentinel->next;

	 if(sentinel->next!=v)
	 perror("Vector not found in bucket...\n");
	 else
	 sentinel->next = sentinel->next->next;
	 */

	// optimistic:
	/*
	 if( b->length > 0 ){

	 int middle = b->length / 2;

	 if( b->pointers[middle]->norm > v->norm ) vPos = 0;
	 else vPos = middle+1;

	 }

	 */

	while (b->pointers[vPos] != v && vPos < b->length) {
		vPos++;
	}

	/*
	 while(vPos < b->length-1){
	 b->pointers[vPos] = b->pointers[vPos+1];
	 vPos++;
	 }
	 b->length--;
	 */

	//        if(vPos >= b->length){
	//          perror("Vector not found in bucket...\n");
	//        return 0;//exit(-1);
	//        }
	/*
	 int i = vPos;
	 while(i<b->length-1){
	 b->pointers[i] = b->pointers[i+1];
	 i++;
	 }
	 */

	b->length--;
	b->pointers[vPos] = b->pointers[b->length];

	// return vPos;
	//}
}

/* Compute the inner product between different vectors v and w */
long ip(struct NodeC *v, struct NodeC *w) {

#ifdef SSE
	register long dot = 0;
	register int i = 0;

	const int loopBound = N - 7;

	__m128i vsum, vecPi, vecCi, vecQCi;

	vsum = _mm_set1_epi16(0);

	int16_t * const pA = v->data;
	int16_t * const pB = w->data;

	for (; i < loopBound; i += 8) {
		vecPi = _mm_load_si128((__m128i *) &(pA)[i]);
		vecCi = _mm_load_si128((__m128i *) &(pB)[i]);
		vecQCi = _mm_madd_epi16(vecPi, vecCi);
		vsum = _mm_add_epi32(vsum, vecQCi);
	}

	vsum = _mm_hadd_epi32(vsum, vsum);
	vsum = _mm_hadd_epi32(vsum, vsum);

	dot += (int32_t) _mm_extract_epi32(vsum, 0);

	for (; i < N; i++)
		dot += pA[i] * pB[i];

	return dot;
#else
	long res = 0;
	for (int i = 0; i < N; i++) {
		res += v->data[i] * w->data[i];
	}
	return res;
#endif
}

long lshash(struct NodeC *v, int t) {

	long res = 0;

	/*   int AbsIPMax = INT_MAX;//100000000;
	 int AbsIPMaxPos = 0;
	 int AbsIPMaxSgn = 0;
	 */
	for (int k = 0; k < K; k++) {
		res <<= 1;
		int IP = v->data[AH[t][k][0]] - v->data[AH[t][k][1]];
		if (IP > 0)
			res++;
	}

	// Merge buckets u and 2^K - u - 1
	if (res >= BUCKETS)
		res = 2 * BUCKETS - res - 1;

	return res;
}

inline int max(int a, int b) {

	if (a > b)
		return a;

	return b;
}

void convertFromStr(NTL::vec_ZZ &value, const std::string &s) {
	std::stringstream ss(s);
	ss >> value;
}

mat_ZZ randomizeMatrix(NTL::mat_ZZ A, NTL::vec_ZZ randoms, int seed, int DIM) {
	ZZ one;
	one = 1;
	mat_ZZ newmat, P;
	newmat.SetDims(DIM, DIM);
	P.SetDims(DIM, DIM);
	int newmati, newmatj, randompicker;

	srand(seed);
	for (int i = 0; i < ((DIM * DIM - DIM) / 2); i++) {
		// cout << "random " << rand() << endl;
		newmati = rand() % (DIM - 1);
		newmatj = rand() % (newmati + 1);
		randompicker = rand() % 5;
		newmat[newmati + 1][newmatj] = randoms[randompicker];
	}
	for (int i = 0; i < DIM; i++) {
		newmat[i][i] = one;
		P[i][i] = one;
	}
	mul(newmat, newmat, A);
	int first, second;
	ZZ temp;
	for (int i = 0; i < DIM; i++) {
		first = rand() % DIM;
		second = rand() % DIM;
		for (int j = 0; j < DIM; j++) {
			temp = P[j][first];
			P[j][first] = P[j][second];
			P[j][second] = temp;
		}
	}
	mul(newmat, P, newmat);
	return newmat;
}

void VecZZToListPoint(const NTL::vec_ZZ &v, struct NodeC *p) {

	p->norm = 0;

	for (int i = 0; i < N; ++i) {

		if (sizeof(MY_DATA) == 2) {

			int int_aux;
			short short_aux;
			int_aux = to_int(v[i]);
			short_aux = (short) int_aux;
			p->data[i] = short_aux;

		} else if (sizeof(MY_DATA) == 4)
			p->data[i] = to_int(v[i]);

		p->norm += p->data[i] * p->data[i];
	}
}

void MatLongtoVect(long **A, struct NodeC *p, int no) {

	p->norm = 0;

	for (int i = 0; i < N; ++i) {

		p->data[i] = (short) A[no][i];

		p->norm += p->data[i] * p->data[i];
	}
}

void MatIntFromMatZZ(long **A, NTL::mat_ZZ B) {

	int row, col, rows = B.NumRows(), cols = B.NumCols();

	for (row = 0; row < rows; ++row) {
		for (col = 0; col < cols; ++col) {
			A[row][col] = to_long(B[row][col]);
		}
	}
}

void MatZZFromMatInt(NTL::mat_ZZ B, long **A) {

	int row, col;

	for (row = 0; row < N; ++row) {
		for (col = 0; col < N; ++col) {
			B[row][col] = A[row][col];
		}
	}
}

void MatDoubleFromMatRR(double **A, NTL::mat_RR B) {

	int row, col, rows = B.NumRows(), cols = B.NumCols();

	for (row = 0; row < rows; ++row) {
		for (col = 0; col < cols; ++col) {
			A[row][col] = to_double(B[row][col]);
		}
	}
}

void insertSortedStack(struct NodeC **pvt_stack, struct NodeC *e,
		long *pvt_stack_size) {

	struct NodeC *sentinel = *pvt_stack;

	/*
	 //no elements on the stack
	 if(!(*pvt_stack)){ *pvt_stack = e;
	 (*pvt_stack_size)++; }

	 //at least one, but new is smaller:
	 else if(e->norm < (*pvt_stack)->norm){
	 e->next = (*pvt_stack);
	 (*pvt_stack) = e;
	 (*pvt_stack_size)++;
	 }
	 */

	/* No elements on the stack or at least one, but new smaller */
	if (*pvt_stack == NULL || (*pvt_stack)->norm >= e->norm) {
		e->next = (*pvt_stack);
		(*pvt_stack) = e;
		(*pvt_stack_size)++;
	} else {
		/* Locate the node before the point of insertion */
		while (sentinel->next && (sentinel->next->norm < e->norm)) {
			sentinel = sentinel->next;
		}
		e->next = sentinel->next;
		sentinel->next = e;
		(*pvt_stack_size)++;
	}

	/*
	 //element to insert is bigger than head
	 else{
	 //there is only one elem in the stack (and it must be smaller because last
	 test failed)
	 if(!((*pvt_stack)->next)){
	 (*pvt_stack)->next=e;
	 (*pvt_stack_size)++;
	 }
	 else{
	 //there are more elements
	 while(sentinel->next){
	 if(sentinel->norm < e->norm &&
	 sentinel->next->norm > e->norm){
	 e->next =
	 sentinel->next;
	 sentinel->next = e;
	 (*pvt_stack_size)++;
	 break;
	 }
	 sentinel = sentinel->next;
	 }
	 if(!sentinel->next){ sentinel->next =
	 e;(*pvt_stack_size)++;}
	 }
	 }

	 */

	return;
}
