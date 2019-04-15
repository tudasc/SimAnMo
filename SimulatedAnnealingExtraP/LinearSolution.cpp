#include "stdafx.h"
#include "LinearSolution.h"
#include <random>

LinearSolution::LinearSolution()
{
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < _len; i++) _coefficients[i] = 0.0;
}

LinearSolution::LinearSolution(MeasurementDB* mdb) {

}

LinearSolution::LinearSolution(double* coefficients)
{
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < _len; i++) _coefficients[i] = coefficients[i];
}

LinearSolution::LinearSolution(const LinearSolution& other) {
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < this->_len; i++)
		this->_coefficients[i] = other._coefficients[i];

	_costs = other._costs;
	setRandomID();
}

LinearSolution & LinearSolution::operator= (const LinearSolution & other) {
	if (_len > 0)
		_coefficients = new double[_len];

	for (int i = 0; i < this->_len; i++)
		this->_coefficients[i] = other._coefficients[i];

	_costs = other._costs;
	setRandomID();
	return *this;
}

LinearSolution::~LinearSolution() {

}

LinearSolution LinearSolution::getNeighborSolution() {
	LinearSolution random_sol = LinearSolution(this->get_coefficients());

	std::random_device seeder;
	std::mt19937 engine(seeder());
	std::uniform_int_distribution<int> dist0_1(0, 1);

	// Decide which coefficient to change c_0 or c_1
	int coeff = dist0_1(engine);

	// Change c_2

	// Change c_3
	return random_sol;
}

double LinearSolution::evaluateModelFunctionAt(double p)
{
	double* c = _coefficients; // just to make access brief
	double exp = c[0] + p * c[1];

	double y = pow(2.0, exp);

#ifdef SOLUTION_DEBUG
	std::cout << "f(" << p << ")" << y << std::endl;
#endif

	return y;
}

double LinearSolution::evaluateConstantTermAt(double p) {
	return evaluateModelFunctionAt(p);
}

void LinearSolution::printModelFunction() {
	double * c = _coefficients;
	std::cout << "(ID: " << this->id << ") \t f(p) = " << c[0] << " + p * " << c[1] << std::endl;
}