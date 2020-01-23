#include "ExponentialPolynomSolution.h"
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

ExponentialPolynomSolution::ExponentialPolynomSolution()
{
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < _len; i++) _coefficients[i] = 0.0;
}

ExponentialPolynomSolution::ExponentialPolynomSolution(MeasurementDB* mdb)
{
	double min_c_2 = 0.1;
	double max_c_2 = 1.0;

	double min_c_3 = 0.0;
	double max_c_3 = Configurator::getInstance().std_exp_range;

	int num_threads = Configurator::getInstance().num_threads;
	int thread_id = omp_get_thread_num();

#ifdef USE_NAG
	ParameterEstimator paramest = ParameterEstimator(mdb);
#else
	EigenParameterEstimator paramest = EigenParameterEstimator(mdb);
#endif
	RSSCostCalculator costcalc = RSSCostCalculator(mdb);

	if (_len > 0)
		_coefficients = new double[_len];

	std::random_device seeder;
	std::mt19937 engine(seeder());
	std::uniform_real_distribution<double> distc2(min_c_2, max_c_2);
	std::uniform_real_distribution<double> distc3(min_c_3, max_c_3);

	double start_vals[5] = { 0, 0, distc2(seeder), distc3(seeder), 0 };
	for (int i = 0; i < _len; i++) _coefficients[i] = start_vals[i];

	paramest.estimateParameters(this);
	costcalc.calculateCost(this);

	ExponentialPolynomSolution act_sol = *this;

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

ExponentialPolynomSolution::ExponentialPolynomSolution(double* coefficients)
{
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < _len; i++) {
		if (i == 2 && coefficients[i] >= 0.05)
			_coefficients[2] = 0.05;
		else
			_coefficients[i] = coefficients[i];
	}
}

ExponentialPolynomSolution::ExponentialPolynomSolution(const ExponentialPolynomSolution& other) {
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < this->_len; i++)
		this->_coefficients[i] = other._coefficients[i];

	_costs = other._costs;
	setRandomID();
}

ExponentialPolynomSolution & ExponentialPolynomSolution::operator= (const ExponentialPolynomSolution & other) {
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < this->_len; i++)
		this->_coefficients[i] = other._coefficients[i];

	_costs = other._costs;
	setRandomID();
	return *this;
}

ExponentialPolynomSolution ExponentialPolynomSolution::getNeighborSolution() {
	ExponentialPolynomSolution random_sol = ExponentialPolynomSolution(this->get_coefficients());

	std::random_device seeder;
	std::mt19937 engine(seeder());
	std::uniform_int_distribution<int> dist2_3(2, 3);
	std::uniform_int_distribution<int> dist20(-Configurator::getInstance().std_exp_range, Configurator::getInstance().std_exp_range);
	std::uniform_int_distribution<int> distc_2_3_change(-300, 300);
	std::uniform_int_distribution<int> dist0or1(0, 1);

	// Decide which coefficient to change c_2, c_3 or c_4
	int coeff = dist2_3(engine);

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
		while (!(new_val > 0.05) || !(new_val < 2.0));

		random_sol.updateAt(2, new_val);
	}

	// Change c_3
	else if (coeff == 3) {
		double val = double(dist20(engine)) / 200.0;

		if (random_sol.getAt(3) + val > 0.00)
			random_sol.updateAt(3, random_sol.getAt(3) + val);
	}


	return random_sol;
}

double ExponentialPolynomSolution::evaluateModelFunctionAt(double p, double scale)
{
	double* c = _coefficients; // just to make access brief
	double exp = c[2] * pow(p, c[3]);

	double y = c[0] + c[1] * pow(2.0, exp);

	if (abs(scale) > 1e-5) {
		y += scale * y;
	}

#ifdef SOLUTION_DEBUG
	std::cout << "f(" << p << ") = " << c[0] << " + " << c[1] << " * 2^"
		<< exp << std::endl;
#endif

	return y;
}

double ExponentialPolynomSolution::evaluateConstantTermAt(double p)
{
	double* c = _coefficients; // just to make access brief
	double exp = c[2] * pow(p, c[3]) /* pow(log10(p), c[4])*/;
	double y = pow(2.0, exp);

#ifdef SOLUTION_DEBUG
	std::cout << "varterm(" << p << ") = " << y << " with exp = " << exp << std::endl;
#endif

	return y;
}

void ExponentialPolynomSolution::printModelFunction() {
	double * c = _coefficients;

	std::cout << "(ID: " << this->id << ") \t f(p) = " << c[0] << " + " << c[1] << " * 2^ ("
		<< c[2] << " * p^" << c[3] << " )" << std::endl;
}

std::string ExponentialPolynomSolution::printModelFunctionLatex(double scale, bool powed) const {
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
			+ str_c2 + " * \\x ^ (" + str_c3 + ")) " + "})";
	}
	else {
		std::string act_func = str_c0 + " + " + str_c1 + " * 2 ^ ("
			+ str_c2 + " * \\x ^ (" + str_c3 + "))";
		func += "(\\x, {" + act_func + "+" + std::to_string(scale) + "*(" + act_func + ")" + "})";
	}

	return func;
}

std::string ExponentialPolynomSolution::printModelFunctionLatexShow() const {
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
		+ str_c2 + " * x ^ {" + str_c3 + "}}";
	return func;
}