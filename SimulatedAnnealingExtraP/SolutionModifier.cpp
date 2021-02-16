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
#include "FactorialSolution.h"

template<class SolutionType, class CostCalulatorType, class ParamEstType>
SolutionModifier<SolutionType, CostCalulatorType, ParamEstType>::SolutionModifier(MeasurementDB* mdb) 
	: _mdb(mdb)
{
#ifdef USE_NAG
	this->param_est = new ParameterEstimator(_mdb);
#else
	this->param_est = new ParamEstType(_mdb);
#endif
	
	this->cost_calc = CostCalulatorType(_mdb);
}

template<class SolutionType, class CostCalulatorType, class ParamEstType>
SolutionType SolutionModifier<SolutionType, CostCalulatorType, ParamEstType>::randomModifySolution(SolutionType * sol)
{

	SolutionType solBack = *sol;
	SolutionType newsol = *sol;
	int count = 0;
	do
	{
		// Backtrack
		if (count > 1000) {
			newsol = solBack;
			count = 0;
			//cout << "CHI";
		}

		count++;

		newsol = newsol.getNeighborSolution();
		int ret = param_est->estimateParameters(&newsol);
		
		if (ret > 0 && ret < 100)
			continue;

		cost_calc.calculateCost(&newsol);

		if (newsol.get_costs() < 1e-10)
		{
			newsol.printModelFunction();
			cout << "ERRR" << endl;
			int stop = 1;
			//cin >> stop;
		}



	} while (std::isnan(newsol.get_costs()) || newsol.get_costs() > Configurator::getInstance().max_cost);

	if (newsol.get_costs() < 1e-10)
	{
		newsol.printModelFunction();
		cout << "ERRR" << endl;
		int stop = 1;
		//cin >> stop;
	}

	return newsol;
}

template<class SolutionType, class CostCalulatorType, class ParamEstType>
SolutionModifier<SolutionType, CostCalulatorType, ParamEstType>::~SolutionModifier() {
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
template class SolutionModifier<FactorialSolution, RSSCostCalculator>;

template class SolutionModifier<Solution, nnrRSSCostCalculator>;
template class SolutionModifier<LinearSolution, nnrRSSCostCalculator>;
template class SolutionModifier<ExtraPSolution, nnrRSSCostCalculator>;
template class SolutionModifier<ExponentialSolution, nnrRSSCostCalculator>;
template class SolutionModifier<ExponentialPolynomSolution, nnrRSSCostCalculator>;
template class SolutionModifier<FactorialSolution, nnrRSSCostCalculator>;