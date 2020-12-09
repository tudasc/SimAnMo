#include "ExtraPSolution.h"
#include <math.h>
#include <iostream>
#ifdef USE_NAG
#include "ParameterEstimator.h"
#else
#include "EigenParameterEstimator.h"
#endif
#include "RSSCostCalculator.h"
#include "Configurator.h"
#include <iomanip>

/****************
	Represents a model of form c_0 + c_1 * p ^ (c_2) * log_2 ^ (c_3) (p)
****************/

blaze_rng::xorshf128 ExtraPSolution::m_rng(time(NULL));

ExtraPSolution::ExtraPSolution()
	: AbstractSolution()
{
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < _len; i++) _coefficients[i] = 0.0;
}

ExtraPSolution::ExtraPSolution(MeasurementDB* mdb)
	: AbstractSolution(mdb) {
	is_wrapped = false;
	double min_c_2 = Configurator::getInstance().min_pol_range;
	double max_c_2 = Configurator::getInstance().max_pol_range;

	double min_c_3 = Configurator::getInstance().min_log_range;
	double max_c_3 = Configurator::getInstance().max_log_range;
#ifdef USE_NAG
	ParameterEstimator paramest = ParameterEstimator(mdb, 0.0);
#else
	EigenParameterEstimator paramest = EigenParameterEstimator(mdb);
#endif
	
	RSSCostCalculator costcalc = RSSCostCalculator(mdb);

	if (_len > 0)
		_coefficients = new double[_len];

	double start_vals[4] = { 0, 0, 0, 0 };
	for (int i = 0; i < _len; i++) _coefficients[i] = start_vals[i];

	paramest.estimateParameters(this);
	costcalc.calculateCost(this);
	this->_costs = std::numeric_limits<double>::max();

	std::random_device seeder;
	std::mt19937 engine(seeder());

	std::uniform_real_distribution<double> distc23_pol(min_c_2, max_c_2);
	std::uniform_real_distribution<double> distc23_log(min_c_3, max_c_3);

	ExtraPSolution act_sol = *this;

	do
	{
		double v1 = distc23_pol(seeder);
		//v1 = 4.3541777161633748;
		act_sol.updateAt(2, v1);
		double v2 = distc23_log(seeder);
		//v2 = 0.76860518759366092;
		act_sol.updateAt(3, v2);
		paramest.estimateParameters(&act_sol);
		costcalc.calculateCost(&act_sol);
		//count++;
	} while ((act_sol.get_costs() > this->get_costs()) || std::isnan(act_sol.get_costs()));

	*this = act_sol;
}

ExtraPSolution::ExtraPSolution(double* coefficients)
	: AbstractSolution (coefficients)
{
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < _len; i++) _coefficients[i] = coefficients[i];
}

ExtraPSolution::ExtraPSolution(const ExtraPSolution& other)
	: AbstractSolution(other)
{
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < this->_len; i++)
		this->_coefficients[i] = other._coefficients[i];

	setRandomID();
}

ExtraPSolution & ExtraPSolution::operator= (const ExtraPSolution & other) {
	AbstractSolution::operator=(other);

	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < this->_len; i++)
		this->_coefficients[i] = other._coefficients[i];

	setRandomID();
	return *this;
}

ExtraPSolution ExtraPSolution::getNeighborSolution() {
	ExtraPSolution random_sol = ExtraPSolution(this->get_coefficients());
	std::random_device seeder;
	std::mt19937 engine(seeder());

	std::uniform_int_distribution<int> dist2_3 (2, 5);
	std::uniform_int_distribution<int> distc_2_3_change (-15, 15);
	
	double pol_diff = abs(Configurator::getInstance().max_pol_range - Configurator::getInstance().min_pol_range);
	double log_diff = abs(Configurator::getInstance().max_log_range - Configurator::getInstance().min_log_range);



	// Decide which coefficient to change c_2 or c_3
	int coeff = dist2_3(m_rng);

	// Change c_2
	if (coeff == 2 || coeff >= 4) {
		double change = 0.0;
		double val = random_sol.getAt(2);
		int cnt = 0;
		do {
			double t = (double)distc_2_3_change(m_rng) / 10;
			change = (t / 100.0) * pol_diff;
		} while (!((val + change) >= Configurator::getInstance().min_pol_range && (val + change <= Configurator::getInstance().max_pol_range)));
		val += change;
		random_sol.updateAt(2, val);
	}

	// Change c_3
	if (coeff >= 3) {
		double change = 0.0;
		double val = random_sol.getAt(3);
		int cnt = 0;
		do {
			//change = ((double)distc_2_3_change(m_rng) / 200.0) * abs(val);
			double t = (double)distc_2_3_change(m_rng) / 10;
			change = (t / 100.0) * log_diff;
			//cout << change << " from " << val << " after " << ++cnt << endl;
		} while (!((val + change) >= Configurator::getInstance().min_log_range && (val + change <= Configurator::getInstance().max_log_range)));
		val += change;
		//cout << "now: " << val << endl;
		random_sol.updateAt(3, val);
	}

	return random_sol;
}

double ExtraPSolution::evaluateModelFunctionAt(double p, double scale) {
	double* c = _coefficients; // just to make access brief
 	double y = c[0] + c[1] * pow(p, c[2]) * pow(log2(p), c[3]);

	if (is_wrapped) {
		return pow(2, y);
	}

	else
		return y;
}

bool ExtraPSolution::isConstantModel() {
	if (_coefficients[1] < abs(10e-10))
		return true;
	return false;
}

double ExtraPSolution::evaluateConstantTermAt(double p)
{
	double* c = _coefficients; // just to make access brief
	double y = pow(p, c[2]) * pow(log2(p), c[3]);
	return y;
}

std::string ExtraPSolution::printModelType() {
	return "PolyomialLogarithmical";
}

std::string ExtraPSolution::getModelFunction() {
	double * c = _coefficients;
	std::stringstream strstr;

	strstr << setprecision(10) << "(ID: " << this->id << ") \t f(p) = " << c[0] << " + " << c[1] << " * p ^ "
		<< c[2] << " * log_2 ^ " << c[3] <<  "(p) " << std::endl;

	return strstr.str();
}

std::string ExtraPSolution::printModelFunctionLatex(double scale, bool powed) const {
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
		func += "(\\x, {" + str_c0 + " + " + str_c1 + " * ( "
			+ "\\x ^ " + str_c2 + ") * log2(\\x) ^" + str_c3
			+ ")});";
	}

	else if (powed) {
		func += "(\\x, { " + str_base + "^(" + str_c0 + " + " + str_c1 + " * ( "
			+ "\\x ^ " + str_c2 + ") * log2(\\x) ^" + str_c3
			+ "))});";
	}

	else {
		std::string  act_func = str_c0 + " + " + str_c1 + " * ( "
			+ "\\x ^ " + str_c2 + ") * log2(\\x) ^" + str_c3
			+ ")";
		func += "(\\x, {" + act_func + "+" + std::to_string(scale) + "*(" + act_func + ")" + "})";
	}
	return func;
}

std::string ExtraPSolution::printModelFunctionLatexShow() const {
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

	func += str_c0 + " + " + str_c1 + " * "
		+ "x ^ {" + str_c2 + "}"; 
		
	if (abs(getAt(3)) > 1e-3) {
		func += " * log2(x) ^ {" + str_c3	+ "}";
	}
	return func;
}

std::string ExtraPSolution::printModelFunctionLatexShow(bool set) const {
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
		+ "x ^ {" + str_c2 + "}";

	if (abs(getAt(3)) > 1e-3) {
		func += " * log2(x) ^ {" + str_c3	+ "}";
	}
	func += "}";
	return func;
}
