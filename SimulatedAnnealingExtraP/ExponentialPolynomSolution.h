#ifndef EXPONENTIALPOLYNOMSOLUTION_H
#define EXPONENTIALPOLYNOMSOLUTION_H

#include "AbstractSolution.h"
#include "xorshf128.tcc"

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
	bool isConstantModel();
	std::string printModelType();
	std::string getModelFunction();

	std::string printModelFunctionLatex(double scale = 0.0, bool powed = false) const;
	std::string printModelFunctionLatexShow() const;

protected:
	const int _len = 5;

private:
	static std::random_device   m_rd;
	static blaze_rng::xorshf128 m_rng;
};

#endif
