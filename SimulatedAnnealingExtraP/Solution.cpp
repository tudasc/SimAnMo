#include "Solution.h"
#include "Configurator.h"
#include <math.h>
#ifdef USE_NAG
#include "ParameterEstimator.h"
#else
#include "EigenParameterEstimator.h"
#endif
#include <iostream>
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
	double max_c_3 = 1e0;

	double min_c_4 = 0.2;
	double max_c_4 = 1.0;

	double min_c_5 = 0;
	double max_c_5 = Configurator::getInstance().max_pol_range;

	int num_threads = Configurator::getInstance().num_threads;
	int thread_id = omp_get_thread_num();

	double split_c4_steps = (abs(min_c_4) + abs(max_c_4)) / num_threads;
	double split_c4_min = thread_id * split_c4_steps;
	double split_c4_max = (thread_id + 1) * split_c4_steps;

#ifdef USE_NAG
	ParameterEstimator paramest = ParameterEstimator(mdb);	
#else
	EigenParameterEstimator paramest = EigenParameterEstimator(mdb);
#endif
	RSSCostCalculator costcalc = RSSCostCalculator(mdb);

	if (_len > 0)
		_coefficients = new double[_len];	

	double start_vals[5] = { 0, 0, 1, 0, 0 };
	for (int i = 0; i < _len; i++) _coefficients[i] = start_vals[i];
	

	paramest.estimateParameters(this);
	costcalc.calculateCost(this);

	std::random_device seeder;
	std::mt19937 engine(seeder());
	std::uniform_real_distribution<double> distc2(min_c_3, max_c_3);
	std::uniform_real_distribution<double> distc3(min_c_4, max_c_4);
	std::uniform_real_distribution<double> distc4(min_c_5, max_c_5);

	Solution act_sol = *this;
	act_sol.updateAt(4, distc4(seeder));
	int count = 0;

	do
	{
		act_sol.updateAt(2, distc2(seeder));
		act_sol.updateAt(3, distc3(seeder));

		paramest.estimateParameters(&act_sol);
		costcalc.calculateCost(&act_sol);
		count++;
	} while ((act_sol.get_costs() > this->get_costs()) || std::isnan(act_sol.get_costs()));

	*this = act_sol;
}

Solution::Solution(double* coefficients) 
{ 
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < _len; i++) _coefficients[i] = coefficients[i];
}

Solution::Solution(const Solution& other) {
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < this->_len; i++)
		this->_coefficients[i] = other._coefficients[i];

	_costs = other._costs;
	setRandomID();
}

Solution & Solution::operator= (const Solution & other) {
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < this->_len; i++)
		this->_coefficients[i] = other._coefficients[i];

	_costs = other._costs;
	setRandomID();
	return *this;
}

Solution Solution::getNeighborSolution() {
	Solution random_sol = Solution(this->get_coefficients());

	std::random_device seeder;
	std::mt19937 engine(seeder());
	std::uniform_int_distribution<int> dist2_4(2, 4);
	std::uniform_int_distribution<int> dist20(-Configurator::getInstance().std_exp_range, Configurator::getInstance().std_exp_range);
	std::uniform_int_distribution<int> distc_2_3_change(-300, 300);
	std::uniform_int_distribution<int> dist0or1(0,1);

	// Decide which coefficient to change c_2, c_3 or c_4
	int coeff = dist2_4(engine);

	/*
	Er muss doch in Solution nicht in Random Sol nach coefficients suchen, oder?
	*/

	// Change c_2
	double new_val = -1000;
	if (coeff == 2) {
		do
		{
			int fac;
			if (dist0or1(engine) == 1)
				fac = 1;
			else
				fac = -1;

			double val = random_sol.getAt(2);
			int log_of_val = int(log10(val));
			double change = fac * pow(10, log_of_val - 1);

			new_val = val + change;
			// Check if we break from 10^x to 10^(x-1) or 10^(x+1)
			if (int(log10(val)) > int(log10(new_val))) {
				change = change * 0.1;
				new_val = val + change;
			}
			else {
				// Everything okay
			}
		}

		// Limit the search space to exclude unrealistic results
		while (!(new_val > 0.1) || !(new_val < 1.0));

		random_sol.updateAt(2, new_val);
	}

	// Change c_3
	else if (coeff == 3) {
		double val = double(dist20(engine)) / 200.0;

		if (random_sol.getAt(3) + val > 0.00)
			random_sol.updateAt(3, random_sol.getAt(3) + val);
	}

	// Change c_4
	else if (coeff == 4) {
		double change = 0.0;
		double val = random_sol.getAt(4);

		do {
			double temp = (double)distc_2_3_change(engine);
			change = temp / 100.0;
		} while (!((val + change) >= 0.0 && (val + change <= Configurator::getInstance().max_pol_range)));
	}

	return random_sol;
}

double Solution::evaluateModelFunctionAt(double p, double scale)
{
	double* c = _coefficients; // just to make access brief
	double exp = c[2] * pow(p, c[3]) /* pow(log10(p), c[4])*/;

	double y = c[0] + c[1] * pow(2.0, exp) * pow(p, c[4]);

#ifdef SOLUTION_DEBUG
	std::cout << "f(" << p << ") = " << c[0] << " + " << c[1] << " * exp" 
		<< y << " with exp = " << exp << std::endl;
#endif

	return y;
}

double Solution::evaluateConstantTermAt(double p)
{
	double* c = _coefficients; // just to make access brief
	double exp = c[2] * pow(p, c[3]) /* pow(log10(p), c[4])*/;
	double y = pow(2.0, exp) * pow(p, c[4]);

#ifdef SOLUTION_DEBUG
	std::cout << "varterm(" << p << ") = " << y << " with exp = " << exp << std::endl;
#endif

	return y;
}

void Solution::printModelFunction() {
	double * c = _coefficients;
	//std::cout << "(ID: " << this->id << ") \t f(p) = " << c[0] << " + " << c[1] << " * 2^ ("
	//	<< c[2] << " * p^" << c[3] << " * log^" << c[4] << "(p) ) * p^" << c[4] << std::endl;

	std::cout << "(ID: " << this->id << ") \t f(p) = " << c[0] << " + " << c[1] << " * 2^ ("
		<< c[2] << " * p^" << c[3] << " ) p^ * " << c[4] << std::endl;
}

std::string Solution::printModelFunctionLatex(double scale) const {
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

	streamObj << getAt(4);
	std::string str_c4 = streamObj.str();
	streamObj.str("");

	std::string func = "";
	if (scale < std::abs(1e-5))
	{
		func += "(\\x, {" + str_c0 + " + " + str_c1 + " * 2 ^ ("
			+ str_c2 + " * \\x ^ (" + str_c3 + ")) * \\x ^" + str_c4 + "})";
	}
	else {
		std::string act_func = str_c0 + " + " + str_c1 + " * 2 ^ ("
			+ str_c2 + " * \\x ^ (" + str_c3 + ")) * \\x ^" + str_c4;
		func += "(\\x, {" + act_func + "+" + std::to_string(scale) + "*(" + act_func + ")" + "})";
	}


	return func;
}

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

	streamObj << getAt(4);
	std::string str_c4 = streamObj.str();
	streamObj.str("");

	std::string func = "";
	func += str_c0 + " + " + str_c1 + " * 2 ^ {"
		+ str_c2 + " * x ^ {" + str_c3 + "}} * x ^ {" + str_c4 + "}";
	return func;
}