#ifndef EXTRAPSOLUTION_H
#define EXTRAPSOLUTION_H

#include "AbstractSolution.h"
#include <string>
#include <sstream>

class ExtraPSolution : public AbstractSolution {
public:
	ExtraPSolution();
	ExtraPSolution(MeasurementDB* mdb);
	ExtraPSolution(double* coefficients);
	ExtraPSolution(const ExtraPSolution& other);
	ExtraPSolution & operator= (const ExtraPSolution & other);
	~ExtraPSolution() {}

	ExtraPSolution getNeighborSolution();
	double evaluateModelFunctionAt(double x);
	double evaluateConstantTermAt(double p);
	void printModelFunction();


	std::string printModelFunctionLatex() const;

	std::string printModelFunctionLatexShow() const;

protected:
	const int _len = 4;
};

#endif
