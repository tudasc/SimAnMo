#include "GeneralParameterEstimator.h"
#include "EigenParameterEstimator.h"
#include "raRSDParameterEstimator.h"
#include "RMSEParameterEstimator.h"
#include "Configurator.h"

GeneralParameterEstimator::GeneralParameterEstimator() {
	this->_inputDB = nullptr;
	this->_esti = nullptr;
}

GeneralParameterEstimator::GeneralParameterEstimator(MeasurementDB* mdb) {
	this->_inputDB = mdb;

	if (Configurator::getInstance().param_est_typ == TYPE_EIGENPARAMETER)
		this->_esti = new EigenParameterEstimator(mdb);
	else if (Configurator::getInstance().param_est_typ == TYPE_ALGLIBPARAMETER)
		this->_esti = new raRSDParameterEstimator(mdb);
	else if (Configurator::getInstance().param_est_typ == TYPE_RMSEPARAMETER)
		this->_esti = new RMSEParameterEstimator(mdb);
	else
		this->_esti = new EigenParameterEstimator(mdb);
}

GeneralParameterEstimator::~GeneralParameterEstimator() {
	delete this->_esti;
}

int GeneralParameterEstimator::estimateParameters(AbstractSolution* sol, double newrelerr) {
	return this->_esti->estimateParameters(sol, newrelerr);
}