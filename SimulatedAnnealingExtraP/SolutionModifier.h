#ifndef SOLUTIONMODIFIER_H
#define SOLUTIONMODIFIER_H

#include "Solution.h"
#include "EigenParameterEstimator.h"
#include "MeasurementDB.h"
#include "RSSCostCalculator.h"

template<class SolutionType, class CostCalulatorType>
class SolutionModifier {
public:
	SolutionModifier() { param_est = NULL; }
	SolutionModifier(MeasurementDB* mdb);
	SolutionModifier(const SolutionModifier& other) {}
	~SolutionModifier();

	SolutionType randomModifySolution(SolutionType* sol);

private:
	MeasurementDB * _mdb;
	EigenParameterEstimator* param_est;
	CostCalulatorType cost_calc;
};

#endif