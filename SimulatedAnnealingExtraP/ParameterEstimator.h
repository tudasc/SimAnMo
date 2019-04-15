#ifndef PARAMETERESTIMATOR_H
#define PARAMETERESTIMATOR_H

#include <nag.h>
#include <nag_stdlib.h>
#include <nage02.h>
#include "MeasurementDB.h"
#include "AbstractSolution.h"

class ParameterEstimator {
public:
	ParameterEstimator(MeasurementDB * _mdb, double relerr=0.0);
	ParameterEstimator();
	ParameterEstimator(const ParameterEstimator& other);
	ParameterEstimator & operator= (const ParameterEstimator & other);
	~ParameterEstimator();
	void estimateParameters(AbstractSolution* sol, double newrelerr = 0.0);
	virtual void block() = 0;

private:
	double _relerr;
	double resmax;
	double tol;
	Integer exit_status;
	Integer	i;
	Integer irank;
	Integer iter;
	Integer m;
	Integer n;
	Integer pda;
	NagError fail;
	Nag_OrderType order;

	double *A;
	double *b;
	double *x;

	MeasurementDB* _mdb;
};

#endif
