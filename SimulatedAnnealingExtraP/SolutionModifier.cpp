#include "stdafx.h"
#include "SolutionModifier.h"
#include <random>
#include <math.h>
#include "Configurator.h"
#include "RSSCostCalculator.h"
#include "LinearSolution.h"
#include "ExtraPSolution.h"

template<class SolutionType, class CostCalulatorType>
SolutionModifier<SolutionType, CostCalulatorType>::SolutionModifier(MeasurementDB* mdb) : _mdb(mdb)
{
	this->param_est = new EigenParameterEstimator(_mdb);
	this->cost_calc = CostCalulatorType(_mdb);
}

template<class SolutionType, class CostCalulatorType>
SolutionType SolutionModifier<SolutionType, CostCalulatorType>::randomModifySolution(SolutionType * sol)
{
	SolutionType newsol = sol->getNeighborSolution();
	param_est->estimateParameters(&newsol);
	cost_calc.calculateCost(&newsol);

	return newsol;
}

template<class SolutionType, class CostCalulatorType>
SolutionModifier<SolutionType, CostCalulatorType>::~SolutionModifier() {
	if (param_est == NULL) {
		delete param_est;
		param_est = NULL;
	}
	
}

template class SolutionModifier<Solution, RSSCostCalculator>;
template class SolutionModifier<LinearSolution, RSSCostCalculator>;
template class SolutionModifier<ExtraPSolution, RSSCostCalculator>;