#include "ExponentialPolynomSolution.h"
#include "Configurator.h"
#include <math.h>
#ifdef USE_NAG
#include "ParameterEstimator.h"
#else
#include "GeneralParameterEstimator.h"
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
blaze_rng::xorshf128 ExponentialPolynomSolution::m_rng((uint32_t)time(NULL));


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
	GeneralParameterEstimator paramest = GeneralParameterEstimator(mdb);
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
	int outercnt = 0;
	int give_up_cnt = 0;

	do
	{
		// Backtrack
		if (count == 100) {
			act_sol.updateAt(3, distc3(seeder));
		}

		// If it seems not possible to find a valid solution in this range, then relax the search space
		if (++outercnt == 5) {
			distc3 = std::uniform_real_distribution<double>(min_c_3, max_c_3);
			act_sol.updateAt(3, distc3(seeder));
			//cout << "Relaxing" << endl;
			outercnt = 0;
		}

		act_sol.updateAt(2, distc2(seeder));
		//act_sol.updateAt(3, distc3(seeder));

		paramest.estimateParameters(&act_sol);
		costcalc.calculateCost(&act_sol);
		count++;

		give_up_cnt++;
		if (give_up_cnt > 10000) {
			*this = act_sol;
			break;
		}

	} while (act_sol.get_costs() == std::numeric_limits<double>::max() || std::isnan(act_sol.get_costs()));

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

	//std::random_device seeder;
	//std::mt19937 engine(seeder());
	//std::mt19937 engine(getRandomID());

	//std::uniform_int_distribution<int> dist2_3(2, 3);
	std::uniform_int_distribution<int> distc_2_3_change(-15, 15);

	const double exp_diff = abs(Configurator::getInstance().max_exp_exp_range - Configurator::getInstance().min_exp_exp_range);
	const double coeff_diff = abs(Configurator::getInstance().max_exp_coeff_range - Configurator::getInstance().min_exp_coeff_range);

	std::uniform_int_distribution<int> dist2or3(2, 3);
	//std::uniform_real_distribution<double> distTrial(1e-4, 3e-1);

	int count = 0;

	// What to change	
	double old_c2 = random_sol.getAt(2);
	double old_c3 = random_sol.getAt(3);

	// Both coefficients in one loop
	double new_val_c2 = 0;
	double new_val_c3 = 0;
	do
	{
		// Backtrack
		if (++count == 15) {
			count = 0;
			random_sol.updateAt(2, old_c2);
			random_sol.updateAt(3, old_c3);
		}

		double t = (double)distc_2_3_change(m_rng) / 100;
		double change = (t / 100.0) * coeff_diff;
		new_val_c2 = random_sol.getAt(2) + change;

		t = (double)distc_2_3_change(m_rng) / 100;
		change = (t / 100.0) * exp_diff;
		new_val_c3 = random_sol.getAt(3) + change;

	} while (	new_val_c2 < Configurator::getInstance().min_exp_coeff_range ||
				new_val_c2 > Configurator::getInstance().max_exp_coeff_range ||
				new_val_c3 < Configurator::getInstance().min_exp_exp_range ||
				new_val_c3 > Configurator::getInstance().max_exp_exp_range
		);

	random_sol.updateAt(2, new_val_c2);
	random_sol.updateAt(3, new_val_c3);




	/*
	int choice = dist2or3(m_rng);
	if (choice == 2 || 1==1)
	{
		double new_val_c2 = -1000;
		// Change c_2
		do
		{
			// Backtrack
			count++;
			if (count == 5) {
				count = 0;
				random_sol.updateAt(2, old_c2);
				cout << "R" << endl;
			}

			double t = (double)distc_2_3_change(m_rng) / 50;
			double change = (t / 100.0) * coeff_diff;
			new_val_c2 = random_sol.getAt(2) + change;

		} while (new_val_c2 < Configurator::getInstance().min_exp_coeff_range ||
			new_val_c2 > Configurator::getInstance().max_exp_coeff_range
			);
		
		random_sol.updateAt(2, new_val_c2);
	}

	if (choice == 1 || 1 == 1)
	{
		double new_val_c3 = -1000;
		// Change c_2
		do
		{
			// Backtrack
			count++;
			if (count == 5) {
				count = 0;
				random_sol.updateAt(2, old_c2);
				random_sol.updateAt(3, old_c3);
				cout << "R" << endl;
			}

			double t = (double)distc_2_3_change(m_rng) / 100;
			double change = (t / 100.0) * exp_diff;
			new_val_c3 = random_sol.getAt(3) + change;

		} while (	new_val_c3 < Configurator::getInstance().min_exp_exp_range ||
					new_val_c3 > Configurator::getInstance().max_exp_exp_range
			);

		random_sol.updateAt(3, new_val_c3);
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

bool ExponentialPolynomSolution::isConstantModel() {
	if(_coefficients[1] < abs(10e-12))
		return true;
	return false;
}

std::string ExponentialPolynomSolution::printModelType() { 
	return "ExponentialPolynomial"; 
}

std::string ExponentialPolynomSolution::getModelFunction() {
	double* c = _coefficients;

	std::stringstream strstr;
	strstr << setprecision(10) << "(ID: " << this->id << ") \t f(p) = " << c[0] << " + " << c[1] << " * 2^ ("
		<< c[2] << " * p^" << c[3] << " )" << std::endl;

	return strstr.str();
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
