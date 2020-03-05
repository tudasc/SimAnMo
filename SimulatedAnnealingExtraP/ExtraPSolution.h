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

	void switchtoLinScale() {
		this->is_wrapped = true;
	}

	void switchtoLocScale() {
		this->is_wrapped = false;
	}

	ExtraPSolution getNeighborSolution();
	double evaluateModelFunctionAt(double x, double scale = 0);
	double evaluateConstantTermAt(double p);
	void printModelFunction();


	std::string printModelFunctionLatex(double scale = 0.0, bool powed = false) const;

	std::string printModelFunctionLatexShow() const;
	std::string printModelFunctionLatexShow(bool set) const;

protected:
	const int _len = 4;
	ExtraPSolution* lin_log_sol;
	bool is_wrapped;
};

#endif
