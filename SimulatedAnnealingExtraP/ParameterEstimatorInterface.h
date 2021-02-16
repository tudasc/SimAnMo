#ifndef PARMATERESTIMATORINTERFACE_H
#define PARMATERESTIMATORINTERFACE_H

#include "AbstractSolution.h"

#define ALLFINE 0
#define ERR_X0_INF 1
#define ERR_X1_INF 2
#define ERR_GENERAL_INVALID 3
#define ERR_BLEICOPTIMIZE_CRASH 4

#define WARN_SMALL_X0 101
#define WARN_SMALL_X1 102

class ParameterEstimatorInterface {
public:
	virtual int estimateParameters(AbstractSolution* sol, double newrelerr = 0.0) = 0;
	ParameterEstimatorInterface() {};
};

#endif
