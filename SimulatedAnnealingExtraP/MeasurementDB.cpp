#include "MeasurementDB.h"
#include "Configurator.h"
#include <iostream>
#include <cmath>

using namespace std;

MeasurementDB::MeasurementDB() {
	this->_no_trainingpoints = 0;
	this->_no_measurements = 0;
}

MeasurementDB* MeasurementDB::cloneToLogVersion(MeasurementDB* inputDB) {
	MeasurementDB* newLogDB = new MeasurementDB();	
	int base = Configurator::getInstance().base_for_lin_log;

	// Keep the training data
	auto itx = inputDB->_xvals.begin();
	auto ity = inputDB->_yvals.begin();

	for (;	itx != inputDB->_xvals.end(); ++itx, ++ity) {
		newLogDB->_xvals.push_back((*itx));
		//newLogDB->_yvals.push_back(log2(*ity));
		newLogDB->_yvals.push_back( (log(*ity) / log(base)) );
		newLogDB->_no_trainingpoints++;
	}

	// Keep additional measurements
	auto itmy = inputDB->_ymeasure.begin();

	for (auto itmx = inputDB->_xmeasure.begin(); itmx != inputDB->_xmeasure.end(); itmx++) {
		newLogDB->_xmeasure.push_back((*itmx));
		newLogDB->_yvals.push_back( (log(*itmy) / log(base)) );
		newLogDB->_no_measurements++;
		itmy++;
	}

	return newLogDB;
}

void MeasurementDB::addTrainingPoint(std::pair<double, double>& measurement) {
	this->_xvals.push_back( measurement.first);
	this->_yvals.push_back(measurement.second);
	this->_no_trainingpoints++;

	if ( (int)_xvals.size() != _no_trainingpoints) {
		std::cerr << "ERROR: MeasurementDB is inconsistent.\n";
		exit(-1);
	}
}

void MeasurementDB::addMeasurementPoint(std::pair<double, double>& measurement)
{
	this->_xmeasure.push_back(measurement.first);
	this->_ymeasure.push_back(measurement.second);
	this->_no_measurements++;

	//cout << "Adding " << measurement.first << " and " << measurement.second << endl;
	if ((int)_xmeasure.size() != _no_measurements) {
		std::cerr << "ERROR: MeasurementDB is inconsistent.\n";
		exit(-1);
	}
}

void MeasurementDB::printMeasurementDB() {
	std::cout << "Printing MeasurementDB:" << std::endl;
	for (int i = 0; i < this->_no_trainingpoints; i++) {
		std::cout << "Entry " << i << "--> x:" << _xvals[i]
			<< " / y:" << _yvals[i] << std::endl;
	}
}

std::pair<double, double> MeasurementDB::getPairAt(int pos) {
	if (pos >= this->_no_trainingpoints) {
		std::cerr << "ERROR: Overflow in MeasurementDB access." << std::endl;
		exit(-2);
	}

	return std::pair<double, double>(_xvals[pos], _yvals[pos]);
}

std::pair<double, double> MeasurementDB::getMeasurePairAt(int pos)
{
	if (pos >= this->_no_measurements) {
		std::cerr << "ERROR: Overflow in MeasurementDB access." << std::endl;
		exit(-3);
	}

	return std::pair<double, double>(_xmeasure[pos], _ymeasure[pos]);
}
