#include "TemperatureInitializer.h"
#include <vector>
#include "Solution.h"
#include <random>
#include <math.h>
#include <limits>

#include "SolutionModifier.h"
#include "RSSCostCalculator.h"
#include "nnrRSSCostCalculator.h"
#include "RMSECostCalculator.h"
#include "R2CostCalculator.h"
#include "ExtraPSolution.h"
#include "ExponentialSolution.h"
#include "ExponentialPolynomSolution.h"

#include "TopRankQueue.h"

template<class SolutionType, class CostCalcType>
double TemperatureInitializer<SolutionType, CostCalcType>::estimateInitialCost(int sample_size, int neighbor_size)
{
	SolutionModifier<SolutionType, CostCalcType> solmod = SolutionModifier<SolutionType, CostCalcType>(_mdb);
	std::vector<std::vector<SolutionType>> solutions;
	solutions.reserve(sample_size);
	solutions.resize(sample_size);

	DarmstadtSC::TopRankQueue<SolutionType> queuei = DarmstadtSC::TopRankQueue<SolutionType>(100);

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

			//queuei.insert(sol_copy);
		}
	}
	//queuei.printQueue();
	//exit(-2);
	
	double min = std::numeric_limits<double>::max();
	double max = 0.0;

	for (int i = 0; i < sample_size; i++) {

		for (int j = 1; j < neighbor_size; j++) {
		
			double diff = abs(solutions[i][0].get_costs() - solutions[i][j].get_costs());

			if (min > diff)
				min = diff;

			if (max < diff) {
				max = diff;
				//cout << diff << endl;
			}
		}
	}

	double T_init = -max / (log(0.9));
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

template class TemperatureInitializer<Solution, RSSCostCalculator>;
template class TemperatureInitializer<ExtraPSolution, RSSCostCalculator>;
template class TemperatureInitializer<ExponentialSolution, RSSCostCalculator>;
template class TemperatureInitializer<ExponentialPolynomSolution, RSSCostCalculator>;

template class TemperatureInitializer<Solution, nnrRSSCostCalculator>;
template class TemperatureInitializer<ExtraPSolution, nnrRSSCostCalculator>;
template class TemperatureInitializer<ExponentialSolution, nnrRSSCostCalculator>;
template class TemperatureInitializer<ExponentialPolynomSolution, nnrRSSCostCalculator>;

template class TemperatureInitializer<Solution, RMSECostCalculator>;
template class TemperatureInitializer<ExtraPSolution, RMSECostCalculator>;
template class TemperatureInitializer<ExponentialSolution, RMSECostCalculator>;
template class TemperatureInitializer<ExponentialPolynomSolution, RMSECostCalculator>;

template class TemperatureInitializer<Solution, R2CostCalculator>;
template class TemperatureInitializer<ExtraPSolution, R2CostCalculator>;
template class TemperatureInitializer<ExponentialSolution, R2CostCalculator>;
template class TemperatureInitializer<ExponentialPolynomSolution, R2CostCalculator>;
