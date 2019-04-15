#include "stdafx.h"
#include "ParameterEstimator.h"
#include <iostream>

ParameterEstimator::ParameterEstimator() {
	A = NULL;
	b = NULL;
	x = NULL;
}

ParameterEstimator::ParameterEstimator(MeasurementDB* mdb, double relerr) : _mdb(mdb), _relerr(relerr) {
	this->order = Nag_RowMajor;
	this->tol = 0.0;
	exit_status = 0;
	INIT_FAIL(fail);

	m = mdb->get_size();
	n = 2; // this must be set to maximum of all solution type
	pda = m + 1;

	// Allocate memory for coefficient matrix and vectors
	/*if (!(A = NAG_ALLOC((n + 3) * (m + 1), double)) ||
		!(b = NAG_ALLOC(m, double)) || !(x = NAG_ALLOC(n, double)))*/
	if (!(A = NAG_ALLOC((m + 2) * (n + 2), double)) ||
		!(b = NAG_ALLOC(m, double)) || !(x = NAG_ALLOC(n + 2, double)))
	{
		std::cout << "Allocation failure in ParameterEstimator" << std::endl;
		exit_status = -1;
		return;
	}

	// Fill values in vector of right-hand side
	for (int i = 0; i < m; i++)
		b[i] = _mdb->getPairAt(i).second;

	// Zero the solution vector
	for (int i = 0; i < n + 2; i++)
		x[i] = 0.0;

	// Only outputs of function
	resmax = 0.0;
	i = 0;
	irank = 0;
	iter = 0;
}

ParameterEstimator::ParameterEstimator(const ParameterEstimator& other) {
	A = other.A;
	b = other.b;
	x = other.x;
}

ParameterEstimator & ParameterEstimator::operator=(const ParameterEstimator & other)
{
	A = other.A;
	b = other.b;
	x = other.x;

	return *this;
}

ParameterEstimator::~ParameterEstimator() {
	if (A != NULL) {
		NAG_FREE(A);
		A = NULL;
	}
		
	if (b != NULL) {
		NAG_FREE(b);
		b = NULL;
	}
		
	if (x != NULL) {
		NAG_FREE(x);
		x = NULL;
	}
		
}

void ParameterEstimator::estimateParameters(AbstractSolution* sol, double newrelerr) {
	if (newrelerr != this->_relerr)
		_relerr = newrelerr;

	exit_status = 0;
	INIT_FAIL(fail);

	// Zero the solution vector
	for (int i = 0; i < n; i++)
		x[i] = 0.0;

	// Fill values in vector of right-hand side
	for (int i = 0; i < m; i++)
		b[i] = _mdb->getPairAt(i).second;

	// Only outputs of function
	resmax = 0.0;
	i = 0;
	irank = 0;
	iter = 0;

	// Set values from solution to A and b
	// Run over all equations
	//std::cout << "A: " << std::endl;
	/*for (int i = 1; i <= m; i++) {
		// Run over all unknowns in the line
		A[(1 - 1) * pda + i - 1] = 1.0;
		A[(2 - 1) * pda + i - 1] = sol->evaluateConstantTermAt(_mdb->getPairAt(i-1).first);
		//std::cout << sol->evaluateConstantTermAt(_mdb->getPairAt(i - 1).first) << std::endl; 
	}*/

#define MA(I, J) A[(I-1)*(n+2) + J - 1]

	for (int I = 1; I <= m; I++) {
		for (int J = 1; J <= 2; J++) {
			if (J == 1)
				MA(I, J) = 1.0;
			else
				MA(I, J) = sol->evaluateConstantTermAt(_mdb->getPairAt(I - 1).first);
		}

		// Run over all unknowns in the line
		//A[(i - 1) * (n + 2) + 1 - 1] = 1.0;
		//A[(i - 1) * (n + 2) + 2 - 1] = sol->evaluateConstantTermAt(_mdb->getPairAt(i - 1).first);
		//std::cout << sol->evaluateConstantTermAt(_mdb->getPairAt(i - 1).first) << std::endl; 
	}

	nag_lone_fit(order, m, A, b, n + 2, tol, x, &resmax, &irank, &iter, &fail);
	//nag_linf_fit(order, m, n, A, b, tol, &_relerr, x, &resmax, &irank, &iter, &fail);

	if (fail.code == 2086) {
		// Update Solution
		for (i = 1; i <= n; ++i) {
			//printf("%10.4f", x[i - 1]);
			sol->updateAt(i - 1, x[i - 1]);
		}
	}

	else if (fail.code != NE_NOERROR && fail.code != 2086) {
		std::cout << "Error from nag_linf_fit (e02gcc).\n" << fail.message 
			<< " (Code:" << fail.code << ")" << std::endl;
		exit_status = 1;
		return;
	}

	else {
		// Update Solution
		for (int k = 1; k <= n; ++k) {
			sol->updateAt(k-1, x[k - 1]);
		}
	}
}