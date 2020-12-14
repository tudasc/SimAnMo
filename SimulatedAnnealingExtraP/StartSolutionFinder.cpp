#include "StartSolutionFinder.h"
#include "RSSCostCalculator.h"
#ifdef USE_NAG
#include "ParameterEstimator.h"
#else
#include "EigenParameterEstimator.h"
#endif

#include <math.h>
#include <iostream>
#include <cmath>
#include "LinearSolution.h"
#include "ExtraPSolution.h"
#include "FactorialSolution.h"
#include "ExponentialSolution.h"
#include "ExponentialPolynomSolution.h"
#include "nnrRSSCostCalculator.h"

template <class SolutionType, class CostCalcType>
StartSolutionFinder<SolutionType, CostCalcType>::StartSolutionFinder(MeasurementDB * mdb) : _mdb(mdb)
{
	min_c_3 = 10e-1;
	max_c_3 = 5e0;
	
	min_c_4 = 0.4;
	max_c_4 = 2.0;
}

template <class SolutionType, class CostCalcType>
double StartSolutionFinder<SolutionType, CostCalcType>::findStartSolution(SolutionType * sol)
{
	return this->findStartSolution(sol, 0, 1);
}

template <class SolutionType, class CostCalcType>
double StartSolutionFinder<SolutionType, CostCalcType>::findStartSolution(SolutionType * sol, int threadno, int nothreads)
{
	*sol = SolutionType(_mdb);
	CostCalcType cost_calc = CostCalcType(_mdb);
	cost_calc.calculateCost(sol);

	return sol->get_costs();
}

// Explicit instantiation of all the types the template will be used with, so that the linker will be able to find them as usual
template class StartSolutionFinder<Solution, RSSCostCalculator>;
template class StartSolutionFinder<LinearSolution, RSSCostCalculator>;
template class StartSolutionFinder<ExtraPSolution, RSSCostCalculator>;
template class StartSolutionFinder<ExponentialSolution, RSSCostCalculator>;
template class StartSolutionFinder<ExponentialPolynomSolution, RSSCostCalculator>;
template class StartSolutionFinder<FactorialSolution, RSSCostCalculator>;

template class StartSolutionFinder<Solution, nnrRSSCostCalculator>;
template class StartSolutionFinder<LinearSolution, nnrRSSCostCalculator>;
template class StartSolutionFinder<ExtraPSolution, nnrRSSCostCalculator>;
template class StartSolutionFinder<ExponentialSolution, nnrRSSCostCalculator>;
template class StartSolutionFinder<ExponentialPolynomSolution, nnrRSSCostCalculator>;
template class StartSolutionFinder<FactorialSolution, nnrRSSCostCalculator>;