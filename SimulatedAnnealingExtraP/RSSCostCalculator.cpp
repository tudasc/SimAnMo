#include "RSSCostCalculator.h"
#include <iostream>

RSSCostCalculator::RSSCostCalculator(MeasurementDB * mdb) : CostCalculator(mdb) {
	this->const_model = calculateConstantModel();
}

double RSSCostCalculator::calculateCost(AbstractSolution* sol) {
	calculateMetrics(sol);
	sol->set_costs( this->RSS );
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