#ifndef EIGENPARAMETERESTIMATOR_H
#define EIGENPARAMETERESTIMATOR_H

#include "MeasurementDB.h"
#include "AbstractSolution.h"

class EigenParameterEstimator {
public:
	EigenParameterEstimator();
	EigenParameterEstimator(MeasurementDB* mdb);
	~EigenParameterEstimator();
	void estimateParameters(AbstractSolution* sol, double newrelerr = 0.0);
private:
	MeasurementDB * _mdb;
};

#endif // !EIGENPARAMETERESTIMATOR_H

