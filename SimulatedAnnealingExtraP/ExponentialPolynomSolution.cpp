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
#include <stdio.h>
#include <time.h>
#include "nnrRSSCostCalculator.h"

#include <iomanip>

ExponentialPolynomSolution::ExponentialPolynomSolution()
{
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < _len; i++) _coefficients[i] = 0.0;
}

ExponentialPolynomSolution::ExponentialPolynomSolution(MeasurementDB* mdb)
{
	double min_c_2 = Configurator::getInstance().min_exp_coeff_range;
	double max_c_2 = Configurator::getInstance().max_exp_coeff_range;

	double min_c_3 = Configurator::getInstance().min_exp_exp_range;
	double max_c_3 = Configurator::getInstance().max_exp_exp_range;

	/*int omp_num_threads = 3;
	int conf_num_threads = Configurator::getInstance().num_threads;

	omp_set_num_threads (5);
#pragma omp parallel num_threads( 5 )
	{
#pragma omp single
		{
			omp_num_threads = omp_get_num_threads();
		}
	}

	if(omp_num_threads != Configurator::getInstance().num_threads) {
		cerr << "Severe configuration error since number of threads is set to " << Configurator::getInstance().num_threads
				<< " while OMP-runtime is configured for " << omp_num_threads << endl;
		exit(2000);
	}*/

	int num_threads = Configurator::getInstance().num_threads;
	int thread_id = omp_get_thread_num();

#ifdef USE_NAG
	ParameterEstimator paramest = ParameterEstimator(mdb);
#else
	EigenParameterEstimator paramest = EigenParameterEstimator(mdb);
#endif
	nnrRSSCostCalculator costcalc = nnrRSSCostCalculator(mdb);

	if (_len > 0)
		_coefficients = new double[_len];

	double c_3_span = max_c_3 - min_c_3;
	double c_3_thread_min = min_c_3 + (c_3_span/(double)num_threads) * (double)thread_id;
	double c_3_thread_max = c_3_thread_min + c_3_span/(double)num_threads;

	//cout << "Thread: " << thread_id << " from " << c_3_thread_min << " to " << c_3_thread_max << endl;

	std::random_device seeder;
	std::mt19937 engine(seeder());
	std::uniform_real_distribution<double> distc2(min_c_2, max_c_2);
	//std::uniform_real_distribution<double> distc3(min_c_3, max_c_3);
	std::uniform_real_distribution<double> distc3(c_3_thread_min, c_3_thread_max);

	double start_vals[5] = { 0, 0, distc2(seeder), distc3(seeder), 0 };
	for (int i = 0; i < _len; i++) _coefficients[i] = start_vals[i];

	paramest.estimateParameters(this);
	costcalc.calculateCost(this);

	ExponentialPolynomSolution act_sol = *this;

	int count = 0;
	do
	{
		// Backtrack
		if (count == 5000) {
			act_sol.updateAt(3, distc3(seeder));
		}

		act_sol.updateAt(2, distc2(seeder));
		//act_sol.updateAt(3, distc3(seeder));

		paramest.estimateParameters(&act_sol);
		costcalc.calculateCost(&act_sol);
		count++;

	} while (/*(act_sol.get_costs() > this->get_costs()) || */std::isnan(act_sol.get_costs()));

	*this = act_sol;
}

ExponentialPolynomSolution::ExponentialPolynomSolution(double* coefficients)
{
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < _len; i++) {
		//if (i == 2 && coefficients[i] >= 0.05)
		//	_coefficients[2] = 0.05;
		//else
			_coefficients[i] = coefficients[i];
	}
}

ExponentialPolynomSolution::ExponentialPolynomSolution(const ExponentialPolynomSolution& other) : AbstractSolution(other) {
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < this->_len; i++)
		this->_coefficients[i] = other._coefficients[i];

	_RSS = other._RSS;
	_nnrRSS = other._nnrRSS;
	_anRSS = other._anRSS;
	_cost_calc_type = other._cost_calc_type;
	_costs = other._costs;
	setRandomID();
}

ExponentialPolynomSolution & ExponentialPolynomSolution::operator= (const ExponentialPolynomSolution & other)  {
	AbstractSolution::operator=(other);
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < this->_len; i++)
		this->_coefficients[i] = other._coefficients[i];

	_costs = other._costs;
	_RSS = other._RSS;
	_nnrRSS = other._nnrRSS;
	_anRSS = other._anRSS;
	_cost_calc_type = other._cost_calc_type;
	setRandomID();
	return *this;
}

ExponentialPolynomSolution ExponentialPolynomSolution::getNeighborSolution() {

	ExponentialPolynomSolution random_sol = ExponentialPolynomSolution(this->get_coefficients());

	std::random_device seeder;
	std::mt19937 engine(seeder());
	//std::mt19937 engine(getRandomID());

	//std::uniform_int_distribution<int> dist2_3(2, 3);
	//std::uniform_int_distribution<int> distc_2_3_change(-300, 300);

	std::uniform_int_distribution<int> dist0or1(0, 1);
	std::uniform_real_distribution<double> distTrial(1e-4, 3e-1);
	double new_val_c1 = -1000;
	double new_val_c2 = -1000;

	int count = 0;

	// Do changes to both constants at once

	// What to change
	int choice = dist0or1(engine);
	double old_c2 = random_sol.getAt(2);
	double old_c3 = random_sol.getAt(3);

	//if (choice == 0)
	{
		// Change c_2
		do
		{
			// Backtrack
			count++;
			if (count == 1000) {
				count = 0;
				random_sol.updateAt(2, old_c2);
				random_sol.updateAt(3, old_c3);
				cout << "R" << endl;
			}

			double sign = 1.0;
			if (dist0or1(engine) == 1) {
				sign = -1.0;
			}
			double perc = sign * distTrial(engine);
			// Make sure, that it is not too small
			//perc += (perc / perc) * 0.001;

			//cout << sign << " / ";

			double change = perc * random_sol.getAt(2);
			new_val_c1 = random_sol.getAt(2) + change;


			sign = 1.0;
			if (dist0or1(engine) == 1) {
				sign = -1.0;
			}

			//cout << sign << endl;

			perc = sign * distTrial(engine);
			new_val_c2 = random_sol.getAt(3) + change;



			//if (abs(perc) < 1e-3) {
				//cerr << "Very small: " << perc << endl;
			//}

			//cout << "Perc: " << perc << endl;

		} while (new_val_c1 < Configurator::getInstance().min_exp_coeff_range ||
			new_val_c1 > Configurator::getInstance().max_exp_coeff_range ||
			new_val_c2 < Configurator::getInstance().min_exp_exp_range ||
			new_val_c2 > Configurator::getInstance().max_exp_exp_range
			);


		if (new_val_c1 < Configurator::getInstance().min_exp_coeff_range) {
			int stop = 1;
		}

		if (new_val_c1 > Configurator::getInstance().max_exp_coeff_range) {
			int stop = 1;
		}

		if (new_val_c2 < Configurator::getInstance().min_exp_exp_range) {
			int stop = 1;
		}

		if (new_val_c2 > Configurator::getInstance().max_exp_exp_range) {
			int stop = 1;
		}


		random_sol.updateAt(2, new_val_c1);
		random_sol.updateAt(3, new_val_c2);
	}

	//else
	/*{
		// Change c_3
		do
		{
			// Backtrack
			count++;
			if (count == 1000) {
				count = 0;
				random_sol.updateAt(3, old_c3);
			}

			double sign = 1.0;
			if (dist0or1(engine) == 1) {
				sign = -1.0;
			}

			double perc = sign * distTrial(engine);

			// Make sure, that it is not too small
			//perc += (perc / perc) * 0.001;
			double change = perc * random_sol.getAt(3);
			new_val = random_sol.getAt(3) + change;

			//cout << "Perc: " << perc << endl;

		} while (!(new_val > Configurator::getInstance().min_exp_exp_range && new_val <= Configurator::getInstance().max_exp_exp_range));
		random_sol.updateAt(3, new_val);
	}*/

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

	std::cout << setprecision(10) << "(ID: " << this->id << ") \t f(p) = " << c[0] << " + " << c[1] << " * 2^ ("
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
	if (std::abs(scale) < 1e-5)
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

	streamObj << round(getAt(0) * 100) / 100; //getAt(0);
	std::string str_c0 = streamObj.str();
	streamObj.str("");

	streamObj << getAt(1); // round(getAt(1) * 100) / 100; //getAt(1);
	std::string str_c1 = streamObj.str();
	streamObj.str("");

	streamObj << round(getAt(2) * 100) / 100; //getAt(2);
	std::string str_c2 = streamObj.str();
	streamObj.str("");

	streamObj << round(getAt(3) * 100) / 100; //getAt(3);
	std::string str_c3 = streamObj.str();
	streamObj.str("");

	std::string func = "";
	func += str_c0 + " + " + str_c1 + " * 2 ^ {"
		+ str_c2 + " * x ^ {" + str_c3 + "}}";
	return func;
}
