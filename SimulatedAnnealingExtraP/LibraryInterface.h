#ifndef LIBRARYINTERFACE_H
#define LIBRARYINTERFACE_H

#include "AbstractSolution.h"
#include "VersionInfo.h"
#include "Configurator.h"
#include <map>

namespace SimAnMo {

	std::string getSimAnMoLibVersion();

	/*struct Costs
	{
		double RSS;
		double anrRSS;
	};*/

class FunctionModel {
public:
	FunctionModel();
	FunctionModel(AbstractSolution* insol);

	// Gives the estimate of the model at x-axis value x
	double evaluateModelFunctionAt(double x);

	// A string containing the type of model 
	// (PolyomialLogarithmical or ExponentialPolynomial)
	std::string getTypeOfModelFunction();

	// A string containing the model in a human readable form 
	std::string getModelFunction();

	// If we have a constant model f(x) = a, then it returns true
	bool isConstant();

	double getRSS();
	double getraRSD();

	// Returns the RSS or the raRSD, depending on the configuration
	// of the annealing process
	double getCosts();

private:
	AbstractSolution* sol;
};

FunctionModel findModel(std::map<double, double>& training_points,
	std::map<double, double>& measurement_points, std::string options);


} // End namespace SimAnMo

#endif // !LIBRARYINTERFACE_H
