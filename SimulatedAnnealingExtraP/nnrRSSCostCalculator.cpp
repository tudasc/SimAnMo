#include "nnrRSSCostCalculator.h"
#include <iostream>

nnrRSSCostCalculator::nnrRSSCostCalculator(MeasurementDB * mdb) : CostCalculator(mdb) {
	this->const_model = calculateConstantModel();
}

double nnrRSSCostCalculator::calculateCost(AbstractSolution* sol) {
	calculateMetrics(sol);
	//if (this->nnrRSS <= 1000000.0)
		sol->set_costs(this->nnrRSS);
	//else
	//	sol->set_costs(std::numeric_limits<double>::max());
	
	return sol->get_costs();
}

double nnrRSSCostCalculator::getConstantModel()
{
	return this->const_model;
}

double nnrRSSCostCalculator::calculateConstantModel()
{
	double mean = 0.0;

	for (int i = 0; i < _mdb->get_size(); i++)
	{
		std::pair<double, double> act_point = _mdb->getPairAt(i);
		mean += act_point.second / (double)_mdb->get_size();
	}

	return mean;
}