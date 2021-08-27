#include "GeneralParameterEstimator.h"
#include "EigenParameterEstimator.h"
#include "raRSDParameterEstimator.h"
#include "RMSEParameterEstimator.h"
#include "R2ParameterEstimator.h"
#include "Configurator.h"

GeneralParameterEstimator::GeneralParameterEstimator() {
	this->_inputDB = nullptr;
	this->_esti = nullptr;
}

GeneralParameterEstimator::GeneralParameterEstimator(MeasurementDB* mdb) {
	this->_inputDB = mdb;

	cout << "Created GeneralParameterEstimator with: " 
		<< Configurator::getInstance().param_est_typ <<  endl;

	if (Configurator::getInstance().param_est_typ == TYPE_EIGENPARAMETER)
		this->_esti = new EigenParameterEstimator(mdb);
	else if (Configurator::getInstance().param_est_typ == TYPE_ALGLIBPARAMETER)
		this->_esti = new raRSDParameterEstimator(mdb);
	else if (Configurator::getInstance().param_est_typ == TYPE_RMSEPARAMETER)
		this->_esti = new RMSEParameterEstimator(mdb);
	else if (Configurator::getInstance().param_est_typ == TYPE_R2PARAMETER)
		this->_esti = new R2ParameterEstimator(mdb);
	else
		this->_esti = new EigenParameterEstimator(mdb);
}

GeneralParameterEstimator::~GeneralParameterEstimator() {
	delete this->_esti;
}

int GeneralParameterEstimator::estimateParameters(AbstractSolution* sol, double newrelerr) {
	return this->_esti->estimateParameters(sol, newrelerr);
}