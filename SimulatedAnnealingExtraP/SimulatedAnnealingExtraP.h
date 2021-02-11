#ifndef SIMULATEDANNEALINGEXTRAP_H
#define SIMULATEDANNEALINGEXTRAP_H

#include "LibraryInterface.h"
#include "nnrRSSCostCalculator.h"
#include "EigenParameterEstimator.h"
#include "SimAnMoTools.h"
#include <map>

SimAnMo::FunctionModel findBestModel(std::map<double, double>& training_points,
	std::map<double, double>& measurement_points, std::string options);

int findAModel(std::string mtype, std::string costcaltype);

/*
*/
template <class T>
T calculateFullModelMetrics(double coeff[5]) {

	// Generate MeasurementDB
	MeasurementDB* db = nullptr;
	if (SimAnMo::readInputMeasurementData(db)) {
		cerr << "Exiting due to measurement data reading error." << endl;
		exit(1001);
	}

	// All available points for training data and cost caluclation!
	db->unifyMeasurementsToTraining();
	db->printMeasurementDB();

	// Set the coefficients for this solutions
	T sol = T(coeff);
	cout << "Checking solution: ";
	sol.printModelFunction();

	// Generate a solution with this new database
	nnrRSSCostCalculator refCostCalc = nnrRSSCostCalculator(db);
	refCostCalc.calculateCost(&sol);
	return sol;
}

#endif
