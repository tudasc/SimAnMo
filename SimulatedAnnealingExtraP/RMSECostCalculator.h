#ifndef RMSECOSTCALCULATOR_H
#define RMSECOSTCALCULATOR_H

#include "CostCalculator.h"
#include "MeasurementDB.h"

class RMSECostCalculator : public CostCalculator {
public:
	RMSECostCalculator() : CostCalculator() {

	}
	RMSECostCalculator(MeasurementDB* mdb);
	double calculateCost(AbstractSolution* sol);

	void set_measurements(MeasurementDB* mdb) { _mdb = mdb; }
	MeasurementDB* get_measurements() { return _mdb; }

	double getConstantModel();
	string getCostTypeString() const {
		return "RMSE";
	}


private:
	double calculateConstantModel();

};

#endif
