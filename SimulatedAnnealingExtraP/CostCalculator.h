#ifndef  COSTCALCULATOR_H
#define COSTCALCULATOR_H

#include "Solution.h"
#include "MeasurementDB.h"
#include <string>

using namespace std;

class CostCalculator {
public:
	CostCalculator() {
		this->RSS = std::numeric_limits<double>::max();
		this->nnrRSS = std::numeric_limits<double>::max();
		this->RMSE = std::numeric_limits<double>::max();
		this->R2 = std::numeric_limits<double>::max();
		this->const_model = std::numeric_limits<double>::max();
		this->_mdb = NULL;
	};
	CostCalculator(MeasurementDB* mdb) : _mdb(mdb) {
		this->RSS = std::numeric_limits<double>::max();
		this->nnrRSS = std::numeric_limits<double>::max();
		this->RMSE = std::numeric_limits<double>::max();
		this->R2 = std::numeric_limits<double>::max();
		this->const_model = std::numeric_limits<double>::max();
	}	

	virtual double calculateCost(AbstractSolution* sol) =0;

	virtual string getCostTypeString() const = 0;

	double getRSS() const {
		return RSS;
	}

	double getnnrRSS() const {
		return nnrRSS;
	}

	double getRMSE() const {
		return RMSE;
	}

protected:
	double calculateMetrics(AbstractSolution* sol);
	double calculateR2(MeasurementDB* _mdb, AbstractSolution* sol, double avg);
	MeasurementDB* _mdb;

	double const_model;

	double RSS;
	double nnrRSS;
	double RMSE;
	double R2;
};

#endif // ! COSTCALCULATOR_H

