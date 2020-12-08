#ifndef EXPONENTIALSOLUTION_H
#define EXPONENTIALSOLUTION_H

#include "AbstractSolution.h"

class ExponentialSolution : public AbstractSolution {
public:
	ExponentialSolution();
	ExponentialSolution(MeasurementDB* mdb);
	ExponentialSolution(double* coefficients);
	ExponentialSolution(const ExponentialSolution& other);
	ExponentialSolution & operator= (const ExponentialSolution & other);
	~ExponentialSolution() {}

	ExponentialSolution getNeighborSolution();

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
