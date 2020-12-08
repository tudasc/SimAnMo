#ifndef EXPONENTIALPOLYNOMSOLUTION_H
#define EXPONENTIALPOLYNOMSOLUTION_H

#include "AbstractSolution.h"

class ExponentialPolynomSolution : public AbstractSolution {
public:
	ExponentialPolynomSolution();
	ExponentialPolynomSolution(MeasurementDB* mdb);
	ExponentialPolynomSolution(double* coefficients);
	ExponentialPolynomSolution(const ExponentialPolynomSolution& other);
	ExponentialPolynomSolution & operator= (const ExponentialPolynomSolution & other);
	~ExponentialPolynomSolution() {}

	ExponentialPolynomSolution getNeighborSolution();

	double evaluateModelFunctionAt(double x, double scale = 0);
	double evaluateConstantTermAt(double x);
	std::string printModelType();
	std::string printModelFunction();

	std::string printModelFunctionLatex(double scale = 0.0, bool powed = false) const;
	std::string printModelFunctionLatexShow() const;

protected:
	const int _len = 5;
};

#endif
