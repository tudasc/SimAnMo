#include "GeneralParameterEstimator.h"
#include "EigenParameterEstimator.h"
#include "raRSDParameterEstimator.h"
#include "Configurator.h"

GeneralParameterEstimator::GeneralParameterEstimator() {

}

GeneralParameterEstimator::GeneralParameterEstimator(MeasurementDB* mdb) {
	this->_inputDB = mdb;

	if (Configurator::getInstance().param_est_typ == TYPE_EIGENPARAMETER)
		this->_esti = new EigenParameterEstimator(mdb);
	else if (Configurator::getInstance().param_est_typ == TYPE_ALGLIBPARAMETER)
		this->_esti = new raRSDParameterEstimator(mdb);
	else
		this->_esti = new EigenParameterEstimator(mdb);
}

GeneralParameterEstimator::~GeneralParameterEstimator() {
	delete this->_esti;
}

void GeneralParameterEstimator::estimateParameters(AbstractSolution* sol, double newrelerr) {
	this->_esti->estimateParameters(sol, newrelerr);
}