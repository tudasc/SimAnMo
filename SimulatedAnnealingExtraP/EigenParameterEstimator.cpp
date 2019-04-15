#include "stdafx.h"
#include "EigenParameterEstimator.h"

#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

EigenParameterEstimator::EigenParameterEstimator()
{
}

EigenParameterEstimator::EigenParameterEstimator(MeasurementDB * mdb) : _mdb(mdb)
{
}

EigenParameterEstimator::~EigenParameterEstimator()
{
}

void EigenParameterEstimator::estimateParameters(AbstractSolution * sol, double newrelerr)
{
	int m = _mdb->get_size();
	int n = 2;

	MatrixXd A = MatrixXd(m, n);
	VectorXd b = VectorXd(m);
	VectorXd x = VectorXd(n);

	// Fill values in vector of right-hand side b
	for (int i = 0; i < m; i++)
		b[i] = _mdb->getPairAt(i).second;

	// Zero the solution vector x
	for (int i = 0; i < n; i++)
		x[i] = 0.0;

	// Fill the coefficient matrix
	for (int row = 0; row < m; row++) {
		A(row, 0) = 1.0;
		A(row, 1) = sol->evaluateConstantTermAt(_mdb->getPairAt(row).first);
	}

	x = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);

	// Update Solution
	for (int i = 1; i <= n; ++i) {
		sol->updateAt(i - 1, x[i - 1]);
	}
}
