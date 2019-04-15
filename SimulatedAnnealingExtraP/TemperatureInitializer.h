#ifndef TEMPERATUREINITIALIZER_H
#define TEMPERATUREINITIALIZER_H

#include "MeasurementDB.h"

template<class SolutionType>
class TemperatureInitializer {
public:
	TemperatureInitializer(MeasurementDB* mdb) :_mdb(mdb) {}
	double estimateInitialCost(int sample_size, int neighbor_size);

private:
	MeasurementDB * _mdb;
};

#endif
