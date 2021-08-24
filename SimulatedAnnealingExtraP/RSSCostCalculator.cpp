#include "RSSCostCalculator.h"
#include <iostream>

RSSCostCalculator::RSSCostCalculator(MeasurementDB * mdb) : CostCalculator(mdb) {
	this->const_model = calculateConstantModel();
}

double RSSCostCalculator::calculateCost(AbstractSolution* sol) {
	calculateMetrics(sol);
	if(Configurator::getInstance().costcalc_type=="nnrrsscostcalculator" ||
		Configurator::getInstance().costcalc_type == "rarsdcost")
		sol->set_costs( this->nnrRSS );
	else if (Configurator::getInstance().costcalc_type == "rmsescostcalculator" ||
		Configurator::getInstance().costcalc_type == "rmsecost")
		sol->set_costs(this->RMSE);
	else
		sol->set_costs(this->RSS);
	return sol->get_costs();
}

double RSSCostCalculator::getConstantModel()
{
	cerr << "Calling non-implemented function RSSCostCalculator::getConstantModel()" << endl;
	exit(10);
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