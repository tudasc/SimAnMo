#include "RSSCostCalculator.h"
#include <iostream>

RSSCostCalculator::RSSCostCalculator(MeasurementDB * mdb) : CostCalculator(mdb) {
	this->const_model = calculateConstantModel();
}

double RSSCostCalculator::calculateCost(AbstractSolution* sol) {	
	sol->set_costs(calculateMetrics(sol));
	return sol->get_costs();
}

double RSSCostCalculator::getConstantModel()
{
	return this->const_model;
}

double RSSCostCalculator::calculateConstantModel()
{
	double mean = 0.0;

	for (int i = 0; i < _mdb->get_size(); i++)
	{
		std::pair<double, double> act_point = _mdb->getPairAt(i);
		mean += act_point.second / (double)_mdb->get_size();
	}

	return mean;
}

double RSSCostCalculator::calculateMetrics(AbstractSolution* sol)
{
	double constrsscost = 0.0;
	double rsscosts = 0.0;
	double nnrrsscost = 0.0;
	double constnnrrsscost = 0.0;
	double averagelen = 0;

	for (int i = 0; i < _mdb->get_size(); i++)
	{
		std::pair<double, double> act_point = _mdb->getPairAt(i);
		double f_i = sol->evaluateModelFunctionAt(act_point.first);

		// RSS
		rsscosts += (act_point.second - f_i) * (act_point.second - f_i);
		constrsscost += (act_point.second - this->const_model) * (act_point.second - this->const_model);

		// nnrRSS
		nnrrsscost += (abs(act_point.second - f_i) / abs(act_point.second));
		constnnrrsscost += (abs(act_point.second - const_model) / abs(act_point.second));

#ifdef RSSCOSTCALCULATORDEBUG
		std::cout << "Costs for " << act_point.first << " = "
			<< act_point.second << " - " << f_i << std::endl;

#endif
	}

	nnrrsscost /= _mdb->get_size();
	constnnrrsscost /= _mdb->get_size();

	// Set costs in relation to costs of constant model
	this->RSS = rsscosts;
	this->nnrRSS = nnrrsscost;

	return rsscosts;
	//return nnrrsscost/* / nnrrsscost constnnrrsscost*/;
}
