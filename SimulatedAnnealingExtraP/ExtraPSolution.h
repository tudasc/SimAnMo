#ifndef EXTRAPSOLUTION_H
#define EXTRAPSOLUTION_H

#include "AbstractSolution.h"
#include <string>
#include <sstream>
#include <random>
#include "xorshf128.tcc"

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
	std::string printModelType();
	std::string printModelFunction();

	void updateAt(int pos, double val) { 
		if(abs(val) > 0.1)
			_coefficients[pos] = val; 
		else
			_coefficients[pos] = 0.0;
	}

	std::string printModelFunctionLatex(double scale = 0.0, bool powed = false) const;

	std::string printModelFunctionLatexShow() const;
	std::string printModelFunctionLatexShow(bool set) const;

protected:
	const int _len = 4;
	ExtraPSolution* lin_log_sol;

private:
    static std::random_device   m_rd;
    static blaze_rng::xorshf128 m_rng;
};



#endif
