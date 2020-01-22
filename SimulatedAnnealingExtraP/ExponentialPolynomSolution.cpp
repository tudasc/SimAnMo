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
		if (i == 2 && coefficients[i] >= 0.1)
			_coefficients[2] = coefficients[i];
		else
			_coefficients[i] = 0.1;
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
	std::uniform_int_distribution<int> dist2_3(2, 3); // To choose coefficient 
	std::uniform_int_distribution<int> dist20(-Configurator::getInstance().std_exp_range, Configurator::getInstance().std_exp_range);
	std::uniform_int_distribution<int> distc_2_3_change(-10000, 10000);

	// Decide which coefficient to change c_2, c_3 or c_4
	int coeff = dist2_3(engine);

	// Change c_2
	if (coeff == 2) {
		double val = 0.0;
		do
		{
			double perc = double(distc_2_3_change(engine)) / 100000.0; // Between |10| and |0.01| % 
			val = perc * random_sol.getAt(2);
			//std::cout << "2: " << val << " to " << random_sol.getAt(3) + val << std::endl;
		}
		while (!((random_sol.getAt(2) + val > 0.10) && (random_sol.getAt(2) + val < 1.0)));
		random_sol.updateAt(2, random_sol.getAt(2) + val);		
		//std::cout << "Fin" << std::endl;
	}

	// Change c_3
	else if (coeff == 3) {
		double val = 0.0;

		do {
			double perc = double(distc_2_3_change(engine)) / 500000.0; // Between |10| and |0.01| % 			
			val = perc * random_sol.getAt(3);
		} 
		while (!((random_sol.getAt(3) + val > 0.00) && (random_sol.getAt(3) + val < Configurator::getInstance().std_exp_range)));
		random_sol.updateAt(3, random_sol.getAt(3) + val);
		//std::cout << "Fin" << std::endl;
	}

	if (random_sol.getAt(2) < 0.1) {
		std::cout << "Warning" << std::endl;
		int stop = 0;
		std::cin >> stop;

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