#ifndef R2COSTCALCUTATOR_H
#define R2COSTCALCUTATOR_H

#include "CostCalculator.h"
#include "MeasurementDB.h"

class R2CostCalculator : public CostCalculator {
public:
	R2CostCalculator() : CostCalculator() {
	}
	R2CostCalculator(MeasurementDB* mdb);
	double calculateCost(AbstractSolution* sol);

	void set_measurements(MeasurementDB* mdb) { _mdb = mdb; }
	MeasurementDB* get_measurements() { return _mdb; }

	double getConstantModel();

	string getCostTypeString() const {
		return "R2";
	}

private:
	double calculateConstantModel();
};

#endif
