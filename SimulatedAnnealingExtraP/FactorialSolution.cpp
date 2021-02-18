#include "FactorialSolution.h"
#include <math.h>
#include <iostream>
#include <random>
#ifdef USE_NAG
#include "ParameterEstimator.h"
#else
#include "GeneralParameterEstimator.h"
#endif
#include "nnrRSSCostCalculator.h"
#include "Configurator.h"
#include <iomanip>

#define LOGVARIANT

blaze_rng::xorshf128 FactorialSolution::m_rng( std::random_device{}() );

/****************
Represents a model of form c_0 + c_1 * ( p!) * log2(p)^c_2
****************/

FactorialSolution::FactorialSolution()
	: AbstractSolution()
{
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < _len; i++) _coefficients[i] = 0.0;
	lin_log_sol = nullptr;
}

FactorialSolution::FactorialSolution(MeasurementDB* mdb)
	: AbstractSolution(mdb) {

	double min_c_2 = Configurator::getInstance().min_pol_range;
	double max_c_2 = Configurator::getInstance().max_pol_range;

	double min_c_3 = Configurator::getInstance().min_log_range;
	double max_c_3 = Configurator::getInstance().max_log_range;
#ifdef USE_NAG
	ParameterEstimator paramest = ParameterEstimator(mdb, 0.0);
#else
	GeneralParameterEstimator paramest = GeneralParameterEstimator(mdb);
#endif
	nnrRSSCostCalculator costcalc = nnrRSSCostCalculator(mdb);

	if (_len > 0)
		_coefficients = new double[_len];

	double start_vals[4] = { 0, 0, 0, 0 };
	for (int i = 0; i < _len; i++) _coefficients[i] = start_vals[i];

	//paramest.estimateParameters(this);
	//costcalc.calculateCost(this);
	this->_costs = std::numeric_limits<double>::max();

	//std::random_device seeder;
	//std::mt19937 engine(seeder());

	std::uniform_real_distribution<double> distc23_pol(min_c_2, max_c_2);
	std::uniform_real_distribution<double> distc23_log(min_c_3, max_c_3);

	FactorialSolution act_sol = *this;

	int max_cost_overflow_cnt = 0;

	do
	{
		act_sol.updateAt(2, distc23_pol(m_rng));
		act_sol.updateAt(3, distc23_log(m_rng));

		//cout << "Hanging" << endl;
		//act_sol.printModelFunction();
		//int stop = 1;
		//cin >> stop;
		int ret = paramest.estimateParameters(&act_sol);
		if (ret > 0 && ret < 100) {
			continue;
		}

		costcalc.calculateCost(&act_sol);

		if (act_sol.get_costs() > Configurator::getInstance().max_cost) {
			max_cost_overflow_cnt++;
		}

		// Relaxing condition if it seems impossible to find a good enough solution
		if (max_cost_overflow_cnt == 500) {
#pragma omp critical
			{
				Configurator::getInstance().max_cost *= 2.0;
				max_cost_overflow_cnt = 0;
			}
		}
		

	} while (act_sol.get_costs() > Configurator::getInstance().max_cost
		|| std::isnan(act_sol.get_costs()));

	*this = act_sol;
}

FactorialSolution::FactorialSolution(double* coefficients)
	: AbstractSolution(coefficients)
{
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < _len; i++) _coefficients[i] = coefficients[i];

	lin_log_sol = nullptr;
}

FactorialSolution::FactorialSolution(const FactorialSolution& other)
	: AbstractSolution(other)
{
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < this->_len; i++)
		this->_coefficients[i] = other._coefficients[i];

	setRandomID();
	lin_log_sol = nullptr;
}

FactorialSolution & FactorialSolution::operator= (const FactorialSolution & other) {
	AbstractSolution::operator=(other);
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < this->_len; i++)
		this->_coefficients[i] = other._coefficients[i];

	lin_log_sol = other.lin_log_sol;
	setRandomID();
	return *this;
}

FactorialSolution FactorialSolution::getNeighborSolution() {
	FactorialSolution random_sol = FactorialSolution(this->get_coefficients());
	//std::random_device seeder;
	//std::mt19937 engine(seeder());
	std::uniform_int_distribution<int> dist2_3(2, 3);
	//std::uniform_int_distribution<int> distc_2_3_change(-300, 300);

	//double pos_pol_diff = Configurator::getInstance().max_pol_range - Configurator::getInstance().min_pol_range;
	//double pos_log_diff = Configurator::getInstance().max_log_range - Configurator::getInstance().min_log_range;

	std::uniform_real_distribution<double> dist_pol_change(0.005, 0.03);
	std::uniform_int_distribution<int> sign_for_change(0, 1);
	//std::uniform_real_distribution<double> dist_log_change(pos_log_diff*(-300), pos_log_diff*(300));

	// Decide which coefficient to change c_2 or c_3
	int coeff = 2;// dist2_3(seeder);

	// Change c_2
	if (coeff == 2) {
		double change_c2 = 0.0;
		double change_c3 = 0.0;
		double val_c2 = random_sol.getAt(2);
		double val_c3 = random_sol.getAt(3);

		int cnt = 0;

		do {
			//double temp = (double)distc_2_3_change(engine);
			double sign_pol = 1.0;
			if (sign_for_change(m_rng/*engine*/) == 0) sign_pol = -1.0;
			double temp = sign_pol * dist_pol_change(m_rng/*engine*/);
			change_c2 = temp * abs(random_sol.getAt(2));

#ifdef LOGVARIANT
			sign_pol = 1.0;
			if (sign_for_change(m_rng/*engine*/) == 0) 
				sign_pol = -1.0;
			temp = sign_pol * dist_pol_change(m_rng/*engine*/);
			change_c3 = temp * abs(random_sol.getAt(3));
#endif

			cnt++;
			if (cnt == 10) {
				change_c2 = 0.0;
				change_c3 = 0.0;
				break;
			}

		} while (!	(
				(val_c2 + change_c2) >= Configurator::getInstance().min_pol_range 
			&&	(val_c2 + change_c2) <= Configurator::getInstance().max_pol_range
#ifdef LOGVARIANT
			&& (val_c3 + change_c3) >= Configurator::getInstance().min_pol_range
			&& (val_c3 + change_c3) <= Configurator::getInstance().max_pol_range
#endif
			)
		);

		val_c2 += change_c2;
		val_c3 += change_c3;

		random_sol.updateAt(2, val_c2);
#ifdef LOGVARIANT
		random_sol.updateAt(3, val_c3);
#endif

	}

	// Change c_3
	/*if (coeff == 3) {
		double change = 0.0;
		double val = random_sol.getAt(3);
		do {
			double temp = (double)dist_pol_change(engine);
			change = temp / 500.0;
		} while (!((val + change) >= Configurator::getInstance().min_log_range && (val + change <= Configurator::getInstance().max_log_range)));
		val += change;
		random_sol.updateAt(3, val);
	}*/

	return random_sol;
}

double FactorialSolution::evaluateModelFunctionAt(double p, double scale) {
	double* c = _coefficients; // just to make access brief
	
	double y = 0.0;
#ifndef	LOGVARIANT
	y = c[0] + c[1] * (double)factorial(p) * pow(p, c[2]);	
#else
	y = c[0] + c[1] * (double)factorial((int)p) * pow(p, c[2]) * pow(log2(p), c[3]);
#endif

	if (is_wrapped) {
		return pow(2, y);
	}

	else
		return y;
}

double FactorialSolution::evaluateConstantTermAt(double p)
{
	double* c = _coefficients; // just to make access brief

	double y = 0.0;
#ifndef	LOGVARIANT
	y = (double)factorial((int)p) * pow(p, c[2]);
#else
	y = (double)factorial((int)p) * pow(p, c[2]) * pow(log2(p), c[3]);
#endif
	return y;
}

// Per definition never constant!
bool FactorialSolution::isConstantModel() {
	return false;
}

std::string FactorialSolution::printModelType() {
	return "Factorial";
}

std::string FactorialSolution::getModelFunction() {
	double * c = _coefficients;
	//std::cout << "(ID: " << this->id << ") \t f(p) = " << c[0] << " * fac(p) " << " + " << c[1]
	//	<< " * fac(p)" << " * p ^ " << c[2] << " * log2(p) ^ " << c[3] << std::endl;

	//std factorial
	//std::cout << "(ID: " << this->id << ") \t f(p) = " << c[0] << " + " << c[1]
	//		<< " * fac(p)" << std::endl;

	std::stringstream strstr;
#ifndef	LOGVARIANT
	strstr << setprecision(10) << "(ID: " << this->id << ") \t f(p) = " << c[0] << " + " << c[1]
		<< " * fac(p)" << "* p^(" << c[2] << ")" << std::endl;
#else
	strstr << setprecision(10) << "(ID: " << this->id << ") \t f(p) = " << c[0] << " + " << c[1]
		<< " * fac(p)" << "* p^(" << c[2] << ") * log2^(" << c[3] << ") (p)" << std::endl;
#endif
	return strstr.str();
}

std::string FactorialSolution::printModelFunctionLatex(double scale, bool powed) const {
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

	std::string str_base = std::to_string(Configurator::getInstance().base_for_lin_log);

	std::string func = "";
	if (scale < std::abs(1e-5) && !powed)
	{
		//func += "(\\x, {" + str_c0 + " + " + str_c1 + " * ( "
		//	+ "log2(\\x) ^ " + str_c2 + ") * factorial(\\x)" + "});";

		//std fac
		//func += "(\\x, {" + str_c0 + " + " + str_c1 + " * "
		//	+ "factorial(\\x)" + "});";
#ifdef	LOGVARIANT
		func += "(\\x, {" + str_c0 + " + " + str_c1 + " * "
			+ "factorial(\\x) * \\x^(" + str_c2 + ")" + "*log2(\\x) ^ " + str_c3  + "});";
#else
		func += "(\\x, {" + str_c0 + " + " + str_c1 + " * "
			+ "factorial(\\x) * \\x^(" + str_c2 + ")" + "});";
#endif
	}

	//To Check
	else if (powed) {
		func += "(\\x, { " + str_base + "^(" + str_c0 + " + " + str_c1 + " * ( "
			+ "log2(\\x) ^ " + str_c2 + ") * factorial(\\x) ^" + ")});";
	}

	//To Check
	else {
		std::string  act_func = str_c0 + " + " + str_c1 + " * ( "
			+ "log(\\x) ^ " + str_c2 + ") * factorial(\\x) ^" + ")";
		func += "(\\x, {" + act_func + "+" + std::to_string(scale) + "*(" + act_func + ")" + "})";
	}
	return func;
}

std::string FactorialSolution::printModelFunctionLatexShow() const {
	std::ostringstream streamObj;

	streamObj << round(getAt(0) * 100) / 100; //getAt(0);
	std::string str_c0 = streamObj.str();
	streamObj.str("");

	streamObj << getAt(1); //round(getAt(1) * 100) / 100; //getAt(1);
	std::string str_c1 = streamObj.str();
	streamObj.str("");

	streamObj << round(getAt(2) * 100) / 100; //getAt(2);
	std::string str_c2 = streamObj.str();
	streamObj.str("");

	streamObj << round(getAt(3) * 100) / 100; //getAt(3);
	std::string str_c3 = streamObj.str();
	streamObj.str("");

	std::string str_base = std::to_string(Configurator::getInstance().base_for_lin_log);

	std::string func = "";

	//func += str_c0 + " + " + str_c1 + " * "
	//	+ "log2(x) ^ {" + str_c2 + "}" + " * factorial(p)";
#ifdef	LOGVARIANT
	func += str_c0 + " + " + str_c1 + " * "
		" factorial(p) * x^{" + str_c2 +"}" + " * log2(x) ^ {" + str_c3 + "}";
#else
	func += str_c0 + " + " + str_c1 + " * "
		" factorial(p) * x^{" + str_c2 + "}";
#endif

	return func;
}

std::string FactorialSolution::printModelFunctionLatexShow(bool set) const {
	cerr << "ERR!!!" << endl;
	int stop = 0;
	cin >> stop;
	std::ostringstream streamObj;

	streamObj << round(getAt(0) * 100) / 100; //getAt(0);
	std::string str_c0 = streamObj.str();
	streamObj.str("");

	streamObj << getAt(1); //round(getAt(1) * 100) / 100; //getAt(1);
	std::string str_c1 = streamObj.str();
	streamObj.str("");

	streamObj << round(getAt(2) * 100) / 100; //getAt(2);
	std::string str_c2 = streamObj.str();
	streamObj.str("");

	streamObj << round(getAt(3) * 100) / 100; //getAt(3);
	std::string str_c3 = streamObj.str();
	streamObj.str("");

	std::string str_base = std::to_string(Configurator::getInstance().base_for_lin_log);

	std::string func = "";

	func += str_base + "^{" + str_c0 + " + " + str_c1 + " * "
		+ "log2(x) ^ {" + str_c2 + "}" + " * log2(x) ^ {" + str_c3 + "}";

	func += " * factorial(p)";
	func += "}";
	return func;
}


uint64_t FactorialSolution::factorial(int inp) {
	uint64_t n = inp;
	uint64_t fac = 1;

	for (uint64_t i = 1; i <= n; ++i)
	{
		fac *= i;
	}

	//cout << "returning " << fac << endl;
	return fac;
}
