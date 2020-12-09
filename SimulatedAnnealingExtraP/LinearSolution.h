#ifndef LINEARSOLUTION_H
#define LINEARSOLUTION_H

#include "AbstractSolution.h"

class LinearSolution : public AbstractSolution  {
public:
	LinearSolution();
	LinearSolution(MeasurementDB* mdb);
	LinearSolution(double* coefficients);
	LinearSolution(const LinearSolution& other);
	LinearSolution & operator= (const LinearSolution & other);
	~LinearSolution();

	LinearSolution getNeighborSolution();

	double evaluateModelFunctionAt(double x, double scale = 0);
	bool isConstantModel();
	double evaluateConstantTermAt(double p);

	std::string printModelType();

	std::string getModelFunction();

	std::string printModelFunctionLatex(double scale=0.0, bool powed = false) const {
		return std::string();
	}

	std::string printModelFunctionLatexShow() const {
		return std::string();
	}

protected:
	const int _len = 2;
};

#endif // !LINEARSOLUTION_H

