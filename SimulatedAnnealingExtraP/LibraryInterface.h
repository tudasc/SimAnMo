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

	double evaluateModelFunctionAt(double x);

	std::string getTypeOfModelFunction();
	std::string getModelFunction();
	bool isConstant();

	double getRSS();
	double getraRSD();
	double getCosts();

private:
	AbstractSolution* sol;
};

FunctionModel findModel(std::map<double, double>& training_points,
	std::map<double, double>& measurement_points, std::string options);


} // End namespace SimAnMo

#endif // !LIBRARYINTERFACE_H
