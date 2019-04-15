#ifndef STARTSOLUTIONFINDER_H
#define STARTSOLUTIONFINDER_H

#include "MeasurementDB.h"
#include "Solution.h"

template <class SolutionType, class CostCalcType>
class StartSolutionFinder {
public:
	StartSolutionFinder() {}
	StartSolutionFinder(MeasurementDB* mdb);

	double findStartSolution(SolutionType* sol);
	double findStartSolution(SolutionType* sol, int threadno, int nothreads);

private:
	MeasurementDB * _mdb;
	double min_c_3, max_c_3;
	double min_c_4, max_c_4;
};

#endif // STARTSOLUTIONFINDER_H

