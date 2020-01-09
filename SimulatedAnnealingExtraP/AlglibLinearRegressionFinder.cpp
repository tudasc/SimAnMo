#include "AlglibLinearRegressionFinder.h"

#include "alglib/dataanalysis.h"
#include <math.h>

using namespace alglib;
LinearSolution AlglibLinearRegressionFinder::findSolution()
{
	int n;

	n = _mdb->get_size();
	ae_int_t info;
	ae_int_t nvars;
	linearmodel model;
	lrreport rep;
	real_1d_array c;

	// Allocate and fill for AlgLib
	real_2d_array xy; 
	xy.setlength(n, 2);

	for (int i = 0; i < n; i++) {
		const std::pair<double, double>& act_pair = _mdb->getPairAt(i);
		xy[i][0] = act_pair.first;
		xy[i][1] = log2(act_pair.second);
	}

	lrbuild(xy, n, 1, info, model, rep);

	if (info != 1) {
		printf("Error from AlgLib Regression Finder (lrbuild).\n");
		return LinearSolution();
	}

	lrunpack(model, c, nvars);

	double coeffs[2] = { c[1], c[0] };
	LinearSolution lin_sol = LinearSolution(coeffs);
	return lin_sol;
}
