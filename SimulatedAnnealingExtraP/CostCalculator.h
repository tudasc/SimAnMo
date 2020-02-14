#ifndef  COSTCALCULATOR_H
#define COSTCALCULATOR_H

#include "Solution.h"
#include "MeasurementDB.h"
#include <string>

using namespace std;

class CostCalculator {
public:
	CostCalculator() {};
	CostCalculator(MeasurementDB* mdb) : _mdb(mdb) {}	

	virtual double calculateCost(AbstractSolution* sol) =0;

	virtual string getCostTypeString() const = 0;

	double getRSS() const {
		return RSS;
	}

	double getnnrRSS() const {
		return nnrRSS;
	}

protected:
	double calculateMetrics(AbstractSolution* sol);
	MeasurementDB* _mdb;

	double const_model;

	double RSS;
	double nnrRSS;
};

#endif // ! COSTCALCULATOR_H

