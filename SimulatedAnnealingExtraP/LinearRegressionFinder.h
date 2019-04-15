#ifndef LINEARREGRESSIONFINDER_H
#define LINEARREGRESSIONFINDER_H

#include "MeasurementDB.h"
#include "LinearSolution.h"

class LinearRegressionFinder {
public:
	LinearRegressionFinder(MeasurementDB* mdb) : _mdb(mdb) {

	}

	LinearSolution findSolution();

private:
	MeasurementDB * _mdb;
};

#endif // !LINEARREGRESSIONFINDER_H

