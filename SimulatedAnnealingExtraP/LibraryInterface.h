#ifndef LIBRARYINTERFACE_H
#define LIBRARYINTERFACE_H

#include "AbstractSolution.h"
#include "Configurator.h"
#include <map>

namespace SimAnMo {

	struct Costs
	{
		double RSS;
		double anrRSS;
	};

class FunctionModel {
public:
	FunctionModel();
	double evaluateModelFunctionAt(double x);

	std::string getTypeOfModelFunction();
	std::string getModelFunction();

	Costs getCosts();


private:
	Costs costs;
	AbstractSolution* sol;
};

FunctionModel findModel(std::map<double, double>& training_points,
	std::map<double, double>& measurement_points, std::string options);



} // End namespace SimAnMo

#endif // !LIBRARYINTERFACE_H