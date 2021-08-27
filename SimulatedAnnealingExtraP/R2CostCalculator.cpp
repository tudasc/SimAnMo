#include "R2CostCalculator.h"
#include <iostream>

R2CostCalculator::R2CostCalculator(MeasurementDB* mdb) : CostCalculator(mdb) {
	this->const_model = calculateConstantModel();
}

double R2CostCalculator::calculateCost(AbstractSolution* sol) {
	calculateMetrics(sol);
	if (Configurator::getInstance().costcalc_type == "r2scostcalculator" ||
		Configurator::getInstance().costcalc_type == "r2cost")
		sol->set_costs(1.0 - this->R2); // We have to maximize R2 not to minimize
	else
		sol->set_costs(this->RSS);
	return sol->get_costs();
}

double R2CostCalculator::getConstantModel()
{
	cerr << "Calling non-implemented function RSSCostCalculator::getConstantModel()" << endl;
	exit(10);
	return this->const_model;
}

double R2CostCalculator::calculateConstantModel()
{
	double mean = 0.0;

	for (int i = 0; i < _mdb->get_size(); i++)
	{
		std::pair<double, double> act_point = _mdb->getPairAt(i);

		mean += act_point.second / (double)_mdb->get_size();
	}

	return mean;
}