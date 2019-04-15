#ifndef RSSCOSTCALCULATOR_H
#define RSSCOSTCALCULATOR_H

#include "CostCalculator.h"
#include "MeasurementDB.h"

class RSSCostCalculator : public CostCalculator {
public:
	RSSCostCalculator() : CostCalculator() {}
	RSSCostCalculator(MeasurementDB * mdb);
	double calculateCost(AbstractSolution* sol);

	void set_measurements(MeasurementDB * mdb) { _mdb = mdb; }
	MeasurementDB * get_measurements() { return _mdb; }

	double getConstantModel();

private:
	double calculateConstantModel();
	double calculateMetrics(AbstractSolution* sol);
	double const_model;

	double RSS;
	double nnrRSS;
};

#endif // !RSSCostCalculator

