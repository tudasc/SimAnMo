#include "stdafx.h"
#include "MeasurementDB.h"
#include <iostream>


MeasurementDB::MeasurementDB() {
	this->_no_trainingpoints = 0;
	this->_no_measurements = 0;
}

void MeasurementDB::addTrainingPoint(std::pair<double, double>& measurement) {
	this->_xvals.push_back( measurement.first);
	this->_yvals.push_back(measurement.second);
	this->_no_trainingpoints++;

	if (_xvals.size() != _no_trainingpoints) {
		std::cerr << "ERROR: MeasurementDB is inconsistent.\n";
		exit(-1);
	}
}

void MeasurementDB::addMeasurementPoint(std::pair<double, double>& measurement)
{
	this->_xmeasure.push_back(measurement.first);
	this->_ymeasure.push_back(measurement.second);
	this->_no_measurements++;

	if (_xmeasure.size() != _no_measurements) {
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
