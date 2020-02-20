#ifndef NNRRSSCOSTCALCULATOR_H
#define NNRRSSCOSTCALCULATOR_H

#include "CostCalculator.h"
#include "MeasurementDB.h"

class nnrRSSCostCalculator : public CostCalculator {
public:
	nnrRSSCostCalculator() : CostCalculator() {

	}
	nnrRSSCostCalculator(MeasurementDB * mdb);
	double calculateCost(AbstractSolution* sol);

	void set_measurements(MeasurementDB * mdb) { _mdb = mdb; }
	MeasurementDB * get_measurements() { return _mdb; }

	double getConstantModel();
	string getCostTypeString() const {
		return "raRSD";
	}


private:
	double calculateConstantModel();	
	
};

#endif // !nnrRSSCostCalculator


