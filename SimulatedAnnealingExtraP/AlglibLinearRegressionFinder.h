#ifndef ALGLIB_REGR_H
#define ALGLIB_REGR_H

#include "MeasurementDB.h"
#include "LinearSolution.h"

class AlglibLinearRegressionFinder {
public:
	AlglibLinearRegressionFinder(MeasurementDB* mdb) : _mdb(mdb) {

	}

	LinearSolution findSolution();

private:
	MeasurementDB * _mdb;
};

#endif
