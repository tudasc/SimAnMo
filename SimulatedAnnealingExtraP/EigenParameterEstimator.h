#ifndef EIGENPARAMETERESTIMATOR_H
#define EIGENPARAMETERESTIMATOR_H

#include "MeasurementDB.h"
#include "AbstractSolution.h"
#include "ParameterEstimatorInterface.h"

class EigenParameterEstimator : public ParameterEstimatorInterface {
public:
	EigenParameterEstimator();
	EigenParameterEstimator(MeasurementDB* mdb);
	~EigenParameterEstimator();
	int estimateParameters(AbstractSolution* sol, double newrelerr=0.0);


private:
	MeasurementDB * _mdb;
};

#endif // !EIGENPARAMETERESTIMATOR_H

