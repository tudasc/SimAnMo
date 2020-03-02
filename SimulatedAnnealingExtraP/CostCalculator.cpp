#include "CostCalculator.h"

double CostCalculator::calculateMetrics(AbstractSolution* sol)
{
	double constrsscost = 0.0;
	double rsscosts = 0.0;
	double nnrrsscost = 0.0;
	double anrss = 0.0;
	double constnnrrsscost = 0.0;
	double averagelen = 0;

	for (int i = 0; i < _mdb->get_size(); i++)
	{
		std::pair<double, double> act_point = _mdb->getPairAt(i);
		double f_i = sol->evaluateModelFunctionAt(act_point.first);
		averagelen += act_point.second;

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

	averagelen /= _mdb->get_size();
	nnrrsscost /= _mdb->get_size();
	constnnrrsscost /= _mdb->get_size();

	// Set costs in relation to costs of constant model
	this->RSS = rsscosts;
	this->nnrRSS = nnrrsscost;

	// Set also to solution
	sol->_RSS = rsscosts;
	sol->_nnrRSS = nnrrsscost;
	// anRSS
	sol->_anRSS = (sqrt(rsscosts) / averagelen) / _mdb->get_size();
	sol->_cost_calc_type = getCostTypeString();

	return rsscosts;
	//
}