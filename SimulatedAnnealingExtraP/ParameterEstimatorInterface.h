#ifndef PARMATERESTIMATORINTERFACE_H
#define PARMATERESTIMATORINTERFACE_H

#include "AbstractSolution.h"

class ParameterEstimatorInterface {
public:
	virtual void estimateParameters(AbstractSolution* sol, double newrelerr = 0.0) = 0;
	ParameterEstimatorInterface() {};
};

#endif
