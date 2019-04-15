#ifndef ABSTRACTSOLUTION_H
#define ABSTRACTSOLUTION_H
#include <random>
#include <iostream>
#include "MeasurementDB.h"
#include <string>

class AbstractSolution {
public:
	AbstractSolution() 
	{
		setRandomID();
	}

	/**
		Creates a randomSolution based on the values in the configurator
	**/
	AbstractSolution(MeasurementDB* mdb)
	{
		setRandomID();
	}

	AbstractSolution(double* coefficients) {
		setRandomID();
	}

	AbstractSolution(const AbstractSolution& other) {
		_costs = other._costs;
		setRandomID();
	}

	AbstractSolution & operator= (const AbstractSolution & other) {
		_costs = other._costs;
		setRandomID();
	}

	~AbstractSolution() {
		delete _coefficients;
	}

	double* get_coefficients() { return _coefficients; }
	void set_coefficients(double* coefficients) { 	for (int i = 0; i < _len; i++) _coefficients[i] = coefficients[i];	}
	double getAt(int pos) const {

		return _coefficients[pos];
	}

	double get_costs() { return _costs; }
	void set_costs(double costs) { _costs = costs; }

	virtual double evaluateModelFunctionAt(double x) = 0;
	virtual double evaluateConstantTermAt(double x) = 0;

	void updateAt(int pos, double val) { _coefficients[pos] = val; }

	virtual void printModelFunction() = 0;
	virtual std::string printModelFunctionLatex() const = 0;
	virtual std::string printModelFunctionLatexShow() const = 0;

protected:
	void setCoefficientLength() { 
		std::cout << "Calling setCoefficientLength of superclass\n";
		//_len = -1; 
	}

	void setRandomID() {
		// Set a random id
		std::random_device seeder;
		std::mt19937 engine(seeder());
		std::uniform_int_distribution<int> dist(1, 999999);
		this->id = dist(engine);
	}

	double * _coefficients;
	double _costs;
	const int _len = -1;
	int id;
};

#endif // !ABSTRACTSOLUTION_H

