#ifndef LINEARREGRESSIONFINDER_H
#define LINEARREGRESSIONFINDER_H

#ifdef USE_NAG
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

#endif
#endif // !LINEARREGRESSIONFINDER_H

