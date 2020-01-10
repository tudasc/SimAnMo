#include "LinearRegressionFinder.h"

#ifdef USE_NAG
#include <nag.h>
#include <nag_stdlib.h>
#include <nagg02.h>
#include <math.h>


LinearSolution LinearRegressionFinder::findSolution()
{
	Integer exit_status = 0, n;
	Nag_SumSquare mean;
	Nag_Boolean weight;
	double a, b, df, err_a, err_b, rsq, rss;
	double *wt = 0, *wtptr, *x = 0, *y = 0;
	NagError fail;

	INIT_FAIL(fail);

	n = _mdb->get_size();
	mean = Nag_AboutMean;
	weight = Nag_FALSE;

	if (!(x = NAG_ALLOC(n, double)) ||
		!(y = NAG_ALLOC(n, double)) || !(wt = NAG_ALLOC(n, double)))
	{
		printf("Nag Linear Regression Allocation failure\n");
		exit_status = -1;
		return LinearSolution();
	}

	wtptr = NULL;

	for (int i = 0; i < n; i++) {
		const std::pair<double, double>& act_pair = _mdb->getPairAt(i);
		x[i] = act_pair.first;
		y[i] = log2(act_pair.second);
	}

	nag_simple_linear_regression(mean, n, x, y, wtptr, &a, &b, &err_a, &err_b,
		&rsq, &rss, &df, &fail);
	if (fail.code != NE_NOERROR) {
		printf("Error from nag_simple_linear_regression (g02cac).\n%s\n",
			fail.message);
		exit_status = 1;
		return LinearSolution();
	}

	double coeffs[2] = { a, b };
	LinearSolution lin_sol = LinearSolution(coeffs);
	return lin_sol;	
}
#endif
