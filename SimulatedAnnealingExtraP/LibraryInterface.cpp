#include "LibraryInterface.h"
#include "SimulatedAnnealingExtraP.h"

#include <iostream>

SimAnMo::FunctionModel::FunctionModel() {
	sol = nullptr;
}

SimAnMo::FunctionModel::FunctionModel(AbstractSolution* insol) {
	sol = insol;
}

double SimAnMo::FunctionModel::evaluateModelFunctionAt(double x) {
	return sol->evaluateModelFunctionAt(x);
}

bool SimAnMo::FunctionModel::isConstant() {
	return sol->isConstantModel();
}


double SimAnMo::FunctionModel::getRSS() {
	return sol->_RSS;
}

double SimAnMo::FunctionModel::getraRSD() {
	return sol->_nnrRSS;
}

double SimAnMo::FunctionModel::getCosts() {
	return sol->get_costs();
}

std::string SimAnMo::FunctionModel::getTypeOfModelFunction() {
	return sol->printModelType();
}

std::string SimAnMo::FunctionModel::getModelFunction () {
	return sol->getModelFunction();
}

SimAnMo::FunctionModel SimAnMo::findModel(std::map<double, double>& training_points,
	std::map<double, double>& measurement_points, std::string options) {
	return findBestModel(training_points, measurement_points, options);
}

std::string SimAnMo::getSimAnMoLibVersion() {
    std::cout << "SimAnMo library in version " << Tutorial_VERSION_MAJOR << "."
              << Tutorial_VERSION_MINOR << std::endl;
	return std::to_string(Tutorial_VERSION_MAJOR) + "." + std::to_string(Tutorial_VERSION_MINOR);
}

