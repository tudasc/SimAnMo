#ifndef FACTORIALSOLUTION_H
#define FACTORIALSOLUTION_H

#include "AbstractSolution.h"
#include <string>
#include <sstream>
#include "xorshf128.tcc"

class FactorialSolution : public AbstractSolution {
public:
	FactorialSolution();
	FactorialSolution(MeasurementDB* mdb);
	FactorialSolution(double* coefficients);
	FactorialSolution(const FactorialSolution& other);
	FactorialSolution & operator= (const FactorialSolution & other);
	~FactorialSolution() {}

	void switchtoLinScale() {
		this->is_wrapped = true;
	}

	void switchtoLocScale() {
		this->is_wrapped = false;
	}

	FactorialSolution getNeighborSolution();
	double evaluateModelFunctionAt(double x, double scale = 0);
	double evaluateConstantTermAt(double p);
	bool isConstantModel();
	std::string printModelType();
	std::string getModelFunction();

	std::string printModelFunctionLatex(double scale = 0.0, bool powed = false) const;

	std::string printModelFunctionLatexShow() const;
	std::string printModelFunctionLatexShow(bool set) const;

protected:
	const int _len = 4;
	FactorialSolution* lin_log_sol;

private:
	uint64_t factorial(int inp);
	static std::random_device   m_rd;
	static blaze_rng::xorshf128 m_rng;
};

#endif