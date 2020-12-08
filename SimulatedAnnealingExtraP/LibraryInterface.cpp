#include "LibraryInterface.h"
#include "SimulatedAnnealingExtraP.h"

SimAnMo::FunctionModel::FunctionModel() {

}

double SimAnMo::FunctionModel::evaluateModelFunctionAt(double x) {
	return sol->evaluateModelFunctionAt(x);
}

SimAnMo::Costs SimAnMo::FunctionModel::getCosts() {
	return this->costs;
}

std::string SimAnMo::FunctionModel::getTypeOfModelFunction() {
	return sol->printModelType();
}

std::string SimAnMo::FunctionModel::getModelFunction () {
	return sol->printModelFunction();
}

SimAnMo::FunctionModel SimAnMo::findModel(std::map<double, double>& training_points,
	std::map<double, double>& measurement_points, std::string options) {
	return findBestModel(training_points, measurement_points, options);
}