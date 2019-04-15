#include "stdafx.h"
#include "TemperatureInitializer.h"
#include <vector>
#include "Solution.h"
#include <random>
#include <math.h>
#include <limits>

#include "SolutionModifier.h"
#include "RSSCostCalculator.h"
#include "EigenParameterEstimator.h"
#include "ExtraPSolution.h"

template<class SolutionType>
double TemperatureInitializer<SolutionType>::estimateInitialCost(int sample_size, int neighbor_size)
{
	SolutionModifier<SolutionType, RSSCostCalculator> solmod = SolutionModifier<SolutionType, RSSCostCalculator>(_mdb);
	std::vector<std::vector<SolutionType>> solutions;
	solutions.reserve(sample_size);
	solutions.resize(sample_size);

	// Generate sample_size random solutions
	for (int i = 0; i < sample_size; i++)
	{
		SolutionType sol = SolutionType(_mdb);
		solutions[i].push_back(sol);

		// Create neighborsize neighbors to the actual Solution
		for (int j = 0; j < neighbor_size; j++) {
			SolutionType sol_copy = sol;
			sol_copy = solmod.randomModifySolution(&sol_copy);
			solutions[i].push_back(sol_copy);
		}
	}
	
	double min = std::numeric_limits<double>::max();
	double max = 0.0;

	for (int i = 0; i < sample_size; i++) {
		for (int j = 1; j < neighbor_size; j++) {
			double diff = abs(solutions[i][0].get_costs() - solutions[i][j].get_costs());

			if (min > diff)
				min = diff;

			if (max < diff)
				max = diff;
		}
	}

	double T_init = -max / (log(0.8));
	return T_init;

	/*SolutionModifier<Solution, RSSCostCalculator> solmod = SolutionModifier<Solution, RSSCostCalculator>(_mdb);
	RSSCostCalculator costcal = RSSCostCalculator(_mdb);
	ParameterEstimator paramest = ParameterEstimator(_mdb);
	std::random_device seeder;
	std::mt19937 engine(seeder());

	std::uniform_real_distribution<double> distcoeff(2e-3, 1e1);
	std::uniform_real_distribution<double> distexpo(0.1, 2.0);

	std::vector<std::vector<Solution>> solutions;
	solutions.reserve(sample_size);
	solutions.resize(sample_size);

	// Generate sample_size random solutions
	for (int i = 0; i < sample_size; i++)
	{
		double val1;
		double val2;
		// Create a new random solution
		do {
			val1 = distcoeff(seeder);
			val2 = distexpo(seeder);
			double coeffs[5] = { 0, 0, val1, val2, 0 };
		} while (val1 * pow(90, val2) > 2000);

		double coeffs[5] = { 0, 0, val1, val2, 0 };
		Solution sol = Solution(coeffs);
		paramest.estimateParameters(&sol);
		costcal.calculateCost(&sol);
		solutions[i].push_back(sol);

		// Create neighborsize neighbors to the actual Solution
		for (int j = 0; j < neighbor_size; j++) {
			Solution sol_copy = sol;
			sol_copy = solmod.randomModifySolution(&sol_copy);
			solutions[i].push_back(sol_copy);
		}
	}

	double min = std::numeric_limits<double>::max();
	double max = 0.0;
	
	for (int i = 0; i < sample_size; i++) {
		for (int j = 1; j < neighbor_size; j++) {
			double diff = abs(solutions[i][0].get_costs() - solutions[i][j].get_costs());

			if (min > diff)
				min = diff;

			if (max < diff)
				max = diff;
		}
	}

	double T_init = -max / (log(0.8));
	return T_init;*/
}

template class TemperatureInitializer<Solution>;
template class TemperatureInitializer<ExtraPSolution>;
