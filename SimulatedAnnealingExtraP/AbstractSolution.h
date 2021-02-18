#ifndef ABSTRACTSOLUTION_H
#define ABSTRACTSOLUTION_H
#include <random>
#include <iostream>
#include "MeasurementDB.h"
#include "Configurator.h"
#include <string>

using namespace std;

class AbstractSolution {
public:
	AbstractSolution() 
	{
		setRandomID();
		is_wrapped = false;

		_costs = std::numeric_limits<double>::max();
		_RSS = std::numeric_limits<double>::max();
		_anRSS = std::numeric_limits<double>::max();
		_nnrRSS = std::numeric_limits<double>::max();
	}

	/**
		Creates a randomSolution based on the values in the configurator
	**/
	AbstractSolution(MeasurementDB* mdb)
	{
		setRandomID();
		is_wrapped = false;
		_costs = std::numeric_limits<double>::max();
		_RSS = std::numeric_limits<double>::max();
		_anRSS = std::numeric_limits<double>::max();
		_nnrRSS = std::numeric_limits<double>::max();
	}

	AbstractSolution(double* coefficients) {
		setRandomID();
		is_wrapped = false;
	}

	AbstractSolution(const AbstractSolution& other) {
		_costs = other._costs;
		_RSS = other._RSS;
		_anRSS = other._anRSS;
		_nnrRSS = other._nnrRSS;
		_cost_calc_type = other._cost_calc_type;
		setRandomID();
		is_wrapped = other.is_wrapped;	
	}

	AbstractSolution & operator= (const AbstractSolution & other) {
		_costs = other._costs;
		_RSS = other._RSS;
		_nnrRSS = other._nnrRSS;
		_anRSS = other._anRSS;
		_cost_calc_type = other._cost_calc_type;
		setRandomID();
		is_wrapped = other.is_wrapped;
		return *this; 
	}

	virtual ~AbstractSolution() {
		delete _coefficients;
	}

	double* get_coefficients() { return _coefficients; }
	void set_coefficients(double* coefficients) { 	for (int i = 0; i < _len; i++) _coefficients[i] = coefficients[i];	}
	double getAt(int pos) const {

		return _coefficients[pos];
	}

	double get_costs() { return _costs; }
	void set_costs(double costs) { _costs = costs; }

	virtual double evaluateModelFunctionAt(double x, double scale=0) = 0;
	virtual double evaluateConstantTermAt(double x) = 0;
	virtual bool isConstantModel() = 0;

	virtual void updateAt(int pos, double val) { _coefficients[pos] = val; }

	virtual std::string printModelFunction() {
		std::string str = getModelFunction();
		std::cout << str;
		return str;
	}
	virtual std::string getModelFunction() = 0;

	virtual std::string printModelFunctionLatex(double scale = 0.0, bool powed = false) const = 0;
	virtual std::string printModelFunctionLatexShow() const = 0;
	virtual std::string printModelType() { return "Abstract"; }

	std::string metricDetailsToStr() {
		std::string str = "";
		str += "RSS: " + std::to_string(this->_RSS) + " / ";
		str += "anRSS: " + std::to_string(this->_anRSS) + " / ";
		str += "arNRS: " + std::to_string(this->_nnrRSS);		
		return str;
	}

	double _RSS;
	double _anRSS;
	double _nnrRSS;
	string _cost_calc_type;

protected:
	void setCoefficientLength() { 
		std::cout << "Calling setCoefficientLength of superclass\n";
		//_len = -1; 
	}

	void setRandomID() {
		// Set a random id
		/*std::random_device seeder;
		std::mt19937 engine(seeder());
		std::uniform_int_distribution<int> dist(1, 999999);
		this->id = dist(engine);*/
		this->id = Configurator::getInstance().glob_id++;
	}

	unsigned long long getRandomID() {
		return this->id;
	}

	double * _coefficients;
	double _costs;
	const int _len = -1;
	unsigned long long id;
	bool is_wrapped;

};

#endif // !ABSTRACTSOLUTION_H

