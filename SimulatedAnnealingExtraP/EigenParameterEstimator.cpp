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
	//cout << "Created EigenParameterEstimator" << endl;
}

EigenParameterEstimator::~EigenParameterEstimator()
{
}

int EigenParameterEstimator::estimateParameters(AbstractSolution * sol, double newrelerr)
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
		//double tmp = sol->evaluateConstantTermAt(_mdb->getPairAt(row).first);
		//cout << tmp << endl;
		A(row, 1) = sol->evaluateConstantTermAt(_mdb->getPairAt(row).first);
	}

	// Try with perhaps more stable version
	/*cout << "The solution using the QR decomposition is:\n"
		<< A.colPivHouseholderQr().solve(b) << endl;*/

	// Old Way
	x = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);

	/*cout << "The solution using the SVD decomposition is:\n"
		<< x << endl;*/

	bool x0_small = false;
	bool x1_small = false;

	// Update Solution
	for (int i = 1; i <= n; ++i) {
		if (abs(x[i - 1]) < 10e-11) {
			if (abs(x[i - 1]) > 0.0) 
				x[i - 1] = 10e-11;
			else
				x[i - 1] = -10e-11;	

			if (i == 1) x0_small = true;
			if (i == 2) x1_small = true;
		}

		if (x[i - 1] != x[i - 1]) {
			if (i == 1) return ERR_X0_INF;
			if (i == 2) return ERR_X1_INF;
		}
			
		sol->updateAt(i - 1, x[i - 1]);
	}

	if (x0_small)
		return WARN_SMALL_X0;

	if (x1_small)
		return WARN_SMALL_X1;

	return 0;
}
