#include "SolutionModifier.h"
#include <random>
#include <math.h>
#include "Configurator.h"
#include "RSSCostCalculator.h"
#include "nnrRSSCostCalculator.h"
#include "LinearSolution.h"
#include "ExtraPSolution.h"
#include "ExponentialSolution.h"
#include "ExponentialPolynomSolution.h"

template<class SolutionType, class CostCalulatorType>
SolutionModifier<SolutionType, CostCalulatorType>::SolutionModifier(MeasurementDB* mdb) : _mdb(mdb)
{
#ifdef USE_NAG
	this->param_est = new ParameterEstimator(_mdb);
#else
	this->param_est = new EigenParameterEstimator(_mdb);
#endif
	
	this->cost_calc = CostCalulatorType(_mdb);
}

template<class SolutionType, class CostCalulatorType>
SolutionType SolutionModifier<SolutionType, CostCalulatorType>::randomModifySolution(SolutionType * sol)
{
	SolutionType newsol;
	do
	{
		newsol = sol->getNeighborSolution();
		param_est->estimateParameters(&newsol);
		cost_calc.calculateCost(&newsol);
	} while (std::isnan(newsol.get_costs()));

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
template class SolutionModifier<ExponentialSolution, RSSCostCalculator>;
template class SolutionModifier<ExponentialPolynomSolution, RSSCostCalculator>;


template class SolutionModifier<Solution, nnrRSSCostCalculator>;
template class SolutionModifier<LinearSolution, nnrRSSCostCalculator>;
template class SolutionModifier<ExtraPSolution, nnrRSSCostCalculator>;
template class SolutionModifier<ExponentialSolution, nnrRSSCostCalculator>;
template class SolutionModifier<ExponentialPolynomSolution, nnrRSSCostCalculator>;