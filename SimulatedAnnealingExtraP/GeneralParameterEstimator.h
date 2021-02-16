#ifndef  GENERALPARAMETERESTIMATOR_H
#define GENERALPARAMETERESTIMATOR_H

#include "MeasurementDB.h"
#include "ParameterEstimatorInterface.h"

class GeneralParameterEstimator : public ParameterEstimatorInterface {
public:
	GeneralParameterEstimator();
	GeneralParameterEstimator(MeasurementDB* inputDB);
	~GeneralParameterEstimator();

	int estimateParameters(AbstractSolution* sol, double newrelerr = 0.0);

private:
	MeasurementDB* _inputDB;
	ParameterEstimatorInterface* _esti;
};

#endif // ! GENERALPARAMETERESTIMATOR_H

