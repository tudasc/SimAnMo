#include "Solution.h"
#include "Configurator.h"
#include <math.h>
#ifdef USE_NAG
#include "ParameterEstimator.h"
#else
#include "EigenParameterEstimator.h"
#endif
#include <iostream>
#include <sstream>
#include <random>
#include <omp.h>
#include "Configurator.h"
#include "RSSCostCalculator.h"
#include <sstream>

Solution::Solution() 
{	
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < _len; i++) _coefficients[i] = 0.0;
}

Solution::Solution(MeasurementDB* mdb)
{
	double min_c_3 = 10e-2;
	double max_c_3 = Configurator::getInstance().std_exp_range;

	double min_c_4 = Configurator::getInstance().min_pol_range;
	double max_c_4 = Configurator::getInstance().max_pol_range;

	//int num_threads = Configurator::getInstance().num_threads;
	//int thread_id = omp_get_thread_num();

#ifdef USE_NAG
	ParameterEstimator paramest = ParameterEstimator(mdb);	
#else
	EigenParameterEstimator paramest = EigenParameterEstimator(mdb);
#endif
	RSSCostCalculator costcalc = RSSCostCalculator(mdb);

	if (_len > 0)
		_coefficients = new double[_len];	

	double start_vals[5] = { 0, 0, 10e-2, 0, 0 };
	for (int i = 0; i < _len; i++) _coefficients[i] = start_vals[i];
	

	paramest.estimateParameters(this);
	costcalc.calculateCost(this);

	std::random_device seeder;
	std::mt19937 engine(seeder());
	std::uniform_real_distribution<double> distc2(min_c_3, max_c_3);
	std::uniform_real_distribution<double> distc3(min_c_4, max_c_4);

	Solution act_sol = *this;
	int count = 0;

	do
	{
		act_sol.updateAt(2, distc2(seeder));
		act_sol.updateAt(3, distc3(seeder));

		std::cout << "Snuff: " << distc3(seeder) << std::endl;

		paramest.estimateParameters(&act_sol);
		costcalc.calculateCost(&act_sol);
		count++;
	} while (std::isnan(act_sol.get_costs()));
	*this = act_sol;
}

/*
*/
Solution::Solution(double* coefficients) 
{ 
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < _len; i++) _coefficients[i] = coefficients[i];
}

/*
*/
Solution::Solution(const Solution& other) {
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < this->_len; i++)
		this->_coefficients[i] = other._coefficients[i];

	_costs = other._costs;
	_RSS = other._RSS;
	_nnrRSS = other._nnrRSS;
	_cost_calc_type = other._cost_calc_type;
	setRandomID();
}

/*
*/
Solution & Solution::operator= (const Solution & other) {
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < this->_len; i++)
		this->_coefficients[i] = other._coefficients[i];

	_costs = other._costs;
	_RSS = other._RSS;
	_nnrRSS = other._nnrRSS;
	_cost_calc_type = other._cost_calc_type;
	setRandomID();
	return *this;
}

/*
*/
Solution Solution::getNeighborSolution() {
	Solution random_sol = Solution(this->get_coefficients());

	std::random_device seeder;
	std::mt19937 engine(seeder());
	std::uniform_int_distribution<int> dist2_3(2, 3);
	//std::uniform_int_distribution<int> dist20(-Configurator::getInstance().std_exp_range, Configurator::getInstance().std_exp_range);
	std::uniform_int_distribution<int> distc_2_3_change(-300, 300);
	std::uniform_real_distribution<double> distTrial(-1, 1);
	double new_val = -1000;

	int valtochange = dist2_3(engine);
	// Do changes to both constants at once
	// Change c_2
	if (valtochange == 2)
	{
		do
		{
			double perc = distTrial(engine) / 350.0;
			double change = perc * random_sol.getAt(2);
			new_val = random_sol.getAt(2) + change;

		} while (!((new_val > 5e-2) && (new_val <= Configurator::getInstance().std_exp_range)));
		random_sol.updateAt(2, new_val);
	}

	// Change c_3
	if (valtochange == 3)
	{
		do
		{
			double perc = distTrial(engine) / 350.0;
			double change = perc * random_sol.getAt(3);
			new_val = random_sol.getAt(3) + change;

		} while (!(new_val > Configurator::getInstance().min_pol_range && new_val < Configurator::getInstance().max_pol_range));
		random_sol.updateAt(3, new_val);
	}

	return random_sol;
}

/*
*/
double Solution::evaluateModelFunctionAt(double p, double scale)
{
	double* c = _coefficients; // just to make access brief
	double exp = c[2] * p;

	double y = c[0] + c[1] * pow(2.0, exp) * pow(p, c[3]);
	return y;
}

bool Solution::isConstantModel() {
	if (_coefficients[1] < abs(10e-10))
		return true;
	return false;
}

/*
*/
double Solution::evaluateConstantTermAt(double p)
{
	double* c = _coefficients; // just to make access brief
	double exp = c[2] * p;
	double y = pow(2.0, exp) * pow(p, c[3]);

	return y;
}

std::string Solution::printModelType() {
	return "Solution";
}

/*
*/
std::string Solution::getModelFunction() {
	double * c = _coefficients;

	std::stringstream strstr;
	strstr << "(ID: " << this->id << ") \t f(p) = " << c[0] << " + " << c[1] << " * 2^ ("
		<< c[2] << " * p" << " ) * p^ " << c[3] << std::endl;

	return strstr.str();
}

/*
*/
std::string Solution::printModelFunctionLatex(double scale, bool powed) const {
	std::ostringstream streamObj;

	streamObj << getAt(0);
	std::string str_c0 = streamObj.str();
	streamObj.str("");

	streamObj << getAt(1);
	std::string str_c1 = streamObj.str();
	streamObj.str("");

	streamObj << getAt(2);
	std::string str_c2 = streamObj.str();
	streamObj.str("");

	streamObj << getAt(3);
	std::string str_c3 = streamObj.str();
	streamObj.str("");

	std::string func = "";
	if (scale < std::abs(1e-5))
	{
		func += "(\\x, {" + str_c0 + " + " + str_c1 + " * 2 ^ ("
			+ str_c2 + " * \\x ) * \\x ^" + str_c3 + "})";
	}
	else {
		std::string act_func = str_c0 + " + " + str_c1 + " * 2 ^ ("
			+ str_c2 + " * \\x ) * \\x ^" + str_c3;
		func += "(\\x, {" + act_func + "+" + std::to_string(scale) + "*(" + act_func + ")" + "})";
	}


	return func;
}

/*
*/
std::string Solution::printModelFunctionLatexShow() const {
	std::ostringstream streamObj;

	streamObj << getAt(0);
	std::string str_c0 = streamObj.str();
	streamObj.str("");

	streamObj << getAt(1);
	std::string str_c1 = streamObj.str();
	streamObj.str("");

	streamObj << getAt(2);
	std::string str_c2 = streamObj.str();
	streamObj.str("");

	streamObj << getAt(3);
	std::string str_c3 = streamObj.str();
	streamObj.str("");

	std::string func = "";
	func += str_c0 + " + " + str_c1 + " * 2 ^ {"
		+ str_c2 + " * x} * x ^ {" + str_c3 + "}";
	return func;
}
