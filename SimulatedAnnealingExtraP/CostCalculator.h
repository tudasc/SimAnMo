#ifndef  COSTCALCULATOR_H
#define COSTCALCULATOR_H

#include "Solution.h"
#include "MeasurementDB.h"

class CostCalculator {
public:
	CostCalculator() {};
	CostCalculator(MeasurementDB* mdb) : _mdb(mdb) {}
	virtual double calculateCost(AbstractSolution* sol) = 0;

protected:
	MeasurementDB* _mdb;
};

#endif // ! COSTCALCULATOR_H

