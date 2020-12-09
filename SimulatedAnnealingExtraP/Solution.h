#ifndef SOLUTION_H
#define SOLUTION_H

#include "AbstractSolution.h"
#include "CostCalculator.h"

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
	bool isConstantModel();
	double evaluateConstantTermAt(double x);
	std::string printModelType();
	std::string getModelFunction();

	std::string printModelFunctionLatex(double scale=0.0, bool powed=false) const;
	std::string printModelFunctionLatexShow() const;

protected:
	const int _len = 5;
	double rss = -1;
	double rnnRSS = -1;
};

#endif
