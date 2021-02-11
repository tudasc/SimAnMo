#ifndef SOLUTIONMODIFIER_H
#define SOLUTIONMODIFIER_H

#include "Solution.h"
#ifdef USE_NAG
#include "ParameterEstimator.h"
#else
#include "GeneralParameterEstimator.h"
#endif

#include "MeasurementDB.h"
#include "RSSCostCalculator.h"

template<class SolutionType, class CostCalulatorType, 
	class ParamEstTyp=GeneralParameterEstimator>
class SolutionModifier {
public:
	SolutionModifier() {
		param_est = NULL;
		_mdb = nullptr;
	}
	SolutionModifier(MeasurementDB* mdb);
	SolutionModifier(const SolutionModifier& other) {
	}
	~SolutionModifier();

	SolutionType randomModifySolution(SolutionType* sol);

private:
	MeasurementDB * _mdb;
#ifdef USE_NAG
	ParameterEstimator* param_est;
#else
	ParameterEstimatorInterface* param_est;
	
#endif
	CostCalulatorType cost_calc;
};

#endif