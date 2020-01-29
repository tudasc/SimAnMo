#include "ExponentialSolution.h"
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

ExponentialSolution::ExponentialSolution()
{
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < _len; i++) _coefficients[i] = 0.0;
}

ExponentialSolution::ExponentialSolution(MeasurementDB* mdb)
{
	double min_c_3 = 10e-2;
	double max_c_3 = 4;

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
	std::uniform_real_distribution<double> distc2(min_c_4, max_c_4);

	ExponentialSolution act_sol = *this;
	int count = 0;
	double mincost = std::numeric_limits<double>::max();

	do
	{
		act_sol.updateAt(2, distc2(seeder));
		paramest.estimateParameters(&act_sol);
		costcalc.calculateCost(&act_sol);
		count++;
	} while ((act_sol.get_costs() > mincost) || std::isnan(act_sol.get_costs()));

	*this = act_sol;
}

ExponentialSolution::ExponentialSolution(double* coefficients)
{
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < _len; i++) _coefficients[i] = coefficients[i];
}

ExponentialSolution::ExponentialSolution(const ExponentialSolution& other) {
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < this->_len; i++)
		this->_coefficients[i] = other._coefficients[i];

	_costs = other._costs;
	setRandomID();
}

ExponentialSolution & ExponentialSolution::operator= (const ExponentialSolution & other) {
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < this->_len; i++)
		this->_coefficients[i] = other._coefficients[i];

	_costs = other._costs;
	setRandomID();
	return *this;
}

ExponentialSolution ExponentialSolution::getNeighborSolution() {
	ExponentialSolution random_sol = ExponentialSolution(this->get_coefficients());

	std::random_device seeder;
	std::mt19937 engine(seeder());
	//std::uniform_int_distribution<int> dist20(-Configurator::getInstance().std_exp_range, Configurator::getInstance().std_exp_range);
	std::uniform_int_distribution<int> dist20(-1.0, 1.0);

	// Change c_3
	double ne_wval = 0;
	do {
		double perc = double(dist20(engine)) / 500.0;
		double val = random_sol.getAt(2) * perc;
		ne_wval = random_sol.getAt(2) + val;
	} while (!((ne_wval > 0.00) && (ne_wval < 2.1)));
	random_sol.updateAt(2, ne_wval);
	return random_sol;
}

double ExponentialSolution::evaluateModelFunctionAt(double p, double scale)
{
	double* c = _coefficients; // just to make access brief
	double y = c[0] + c[1] * pow(2.0, c[2]*p);

	if (abs(scale) > 1e-5) {
		y += scale * y;
	}

#ifdef SOLUTION_DEBUG
	std::cout << "f(" << p << ") = " << c[0] << " + " << c[1] << " * exp"
		<< "(" << c[2] << "*p)";
		<< y << std::endl;
#endif

	return y;
}

double ExponentialSolution::evaluateConstantTermAt(double p)
{
	double* c = _coefficients; // just to make access brief
	double y = pow(2.0, c[2] * p);

#ifdef SOLUTION_DEBUG
	std::cout << "varterm(" << p << ") = " << y << std::endl;
#endif

	return y;
}

void ExponentialSolution::printModelFunction() {
	double * c = _coefficients;

	std::cout << "(ID: " << this->id << ") \t f(p) = " << c[0] << " + " << c[1] << " * 2^ ("
		<< c[2] << " * p)" << std::endl;
}

std::string ExponentialSolution::printModelFunctionLatex(double scale, bool powed) const {
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
	if (std::abs(scale) < std::abs(1e-5))
	{
		func += "(\\x, {" + str_c0 + " + " + str_c1 + " * 2 ^ (" + str_c2 + " * \\x )})";
	}

	else
	{
		std::string  act_func = str_c0 + " + " + str_c1 + " * 2 ^ ("+ str_c2 + " * \\x )";
		func += "(\\x, {" + act_func + "+" + std::to_string(scale) + "*(" + act_func +")" + "})";		
	}

	return func;
}

std::string ExponentialSolution::printModelFunctionLatexShow() const {
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
		+ str_c2 + " * x }";
	return func;
}