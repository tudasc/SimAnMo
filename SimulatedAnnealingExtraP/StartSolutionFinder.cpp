#include "stdafx.h"
#include "StartSolutionFinder.h"
#include "RSSCostCalculator.h"
#include "EigenParameterEstimator.h"
#include <math.h>
#include <iostream>
#include <cmath>
#include "LinearSolution.h"
#include "ExtraPSolution.h"

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


	/*CostCalcType costcalc = CostCalcType(_mdb);
	ParameterEstimator paramest = ParameterEstimator(_mdb);

	double start_vals[5] = { 0, 0, 1, 0, 0 };
	sol->set_coefficients(start_vals);
	paramest.estimateParameters(sol);
	costcalc.calculateCost(sol);

	double split_c4_steps = (abs(min_c_4) + abs(max_c_4)) / nothreads;

	double split_c4_min = threadno * split_c4_steps;
	double split_c4_max = (threadno + 1) * split_c4_steps;

	//std::cout << split_c4_min << " / " << split_c4_max << std::endl;

	SolutionType act_sol = *sol;
	bool initial_change = false;

	do {		
		for (double c_4 = split_c4_min; c_4 < split_c4_max; c_4 += 0.1) {
			act_sol.updateAt(3, c_4);
			for (int exp = (int)log10(min_c_3); exp < (int)log10(max_c_3); exp++) {
				int var = exp - 1;
				double change = pow(10, var);
				double c_3 = pow(10, exp);
				for (int i = 0; i < 10; i++) {
					c_3 += change;
					act_sol.updateAt(2, c_3);

					paramest.estimateParameters(&act_sol);
					costcalc.calculateCost(&act_sol);

					if (act_sol.get_costs() < sol->get_costs()) {
						*sol = act_sol;
						initial_change = true;
					}
				}
			}
			split_c4_min = split_c4_min - 0.1;
		} 
	} while (!initial_change);*/

	return sol->get_costs();
}

// Explicit instantiation of all the types the template will be used with, so that the linker will be able to find them as usual
template class StartSolutionFinder<Solution, RSSCostCalculator>;
template class StartSolutionFinder<LinearSolution, RSSCostCalculator>;
template class StartSolutionFinder<ExtraPSolution, RSSCostCalculator>;