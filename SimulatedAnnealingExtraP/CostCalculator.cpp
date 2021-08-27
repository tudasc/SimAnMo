#include "CostCalculator.h"

double CostCalculator::calculateR2(MeasurementDB* _mdb, AbstractSolution* sol, double avg) {
	
	double xtimesy = 0.0;
	double xx = 0.0;
	double yy = 0.0;
	double sx = 0.0;
	double sy = 0.0;

	const double n = (double)_mdb->get_size();

	for (int i = 0; i < n; i++) {
		std::pair<double, double> act_point = _mdb->getPairAt(i);

		sx += act_point.second;
		sy += sol->evaluateModelFunctionAt(act_point.first);
		xx += act_point.second * act_point.second;
		yy += sol->evaluateModelFunctionAt(act_point.first) * sol->evaluateModelFunctionAt(act_point.first);
		xtimesy += act_point.second * sol->evaluateModelFunctionAt(act_point.first);
	}

	double nSxy = n * xtimesy;
	double nominator = nSxy - sx * sy;
	double denominator = sqrt((n * xx - sx * sx) * (n * yy - sy * sy));
	double r = nominator / denominator;
	return r * r;

	double SStot = 0.0;
	double SSres = 0.0;
	double residual = 0.0;
	double residual2 = 0.0;

	for (int i = 0; i < _mdb->get_size(); i++) {
		std::pair<double, double> act_point = _mdb->getPairAt(i);

		// sum up the square of the residual 
		residual = act_point.second - avg;
		SStot = SStot + (residual * residual);

		residual2 = act_point.second - sol->evaluateModelFunctionAt(act_point.first);
		SSres = SSres + (residual2 * residual2);
	}

	return (1.0 - SSres/SStot);
}

double CostCalculator::calculateMetrics(AbstractSolution* sol)
{
	double constrsscost = 0.0;
	double rsscosts = 0.0;
	double nnrrsscost = 0.0;
	double rmsecost = 0.0;
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

		// RMSE
		rmsecost += (act_point.second - f_i) * (act_point.second - f_i);

#ifdef RSSCOSTCALCULATORDEBUG
		std::cout << "Costs for " << act_point.first << " = "
			<< act_point.second << " - " << f_i << std::endl;

#endif
	}

	averagelen /= _mdb->get_size();
	nnrrsscost /= _mdb->get_size();
	constnnrrsscost /= _mdb->get_size();
	rmsecost = sqrt(rmsecost / (double)_mdb->get_size());

	// Set costs in relation to costs of constant model
	this->RSS = rsscosts;
	this->nnrRSS = nnrrsscost;
	this->RMSE = rmsecost;
	this->R2 = calculateR2(_mdb, sol, averagelen);

	// Set also to solution
	sol->_RSS = rsscosts;
	sol->_nnrRSS = nnrrsscost;
	sol->_RMSE = rmsecost;
	// anRSS
	sol->_anRSS = (sqrt(rsscosts) / averagelen) / _mdb->get_size();
	sol->_cost_calc_type = getCostTypeString();

	sol->_R2 = this->R2;

	return rsscosts;
	//
}
