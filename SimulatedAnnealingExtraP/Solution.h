#ifndef SOLUTION_H
#define SOLUTION_H

#include "AbstractSolution.h"

class Solution : public AbstractSolution {
public:
	Solution();
	Solution( MeasurementDB* mdb );
	Solution(double* coefficients);
	Solution(const Solution& other);
	Solution & operator= (const Solution & other);
	~Solution() {}

	Solution getNeighborSolution();

	double evaluateModelFunctionAt(double x, double scale = 0);
	double evaluateConstantTermAt(double x);
	void printModelFunction();	

	std::string printModelFunctionLatex(double scale=0.0) const;
	std::string printModelFunctionLatexShow() const;

protected:
	const int _len = 5;
};

#endif
