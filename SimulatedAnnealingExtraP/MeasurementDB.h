#ifndef MEASUREMENTDB_H
#define MEASUREMENTDB_H

#include <vector>
#include <utility>

class MeasurementDB {
public:
	MeasurementDB();
	void addTrainingPoint(std::pair<double, double>& measurement);
	void addMeasurementPoint(std::pair<double, double>& measurement);

	int get_size() { return this->_no_trainingpoints; }
	int get_measures_size() { return this->_no_measurements; }

	std::pair<double, double> getPairAt(int pos);
	std::pair<double, double> getMeasurePairAt(int pos);

	void printMeasurementDB();

private:
	// Keep the training data
	std::vector<double> _xvals;
	std::vector<double> _yvals;
	int _no_trainingpoints;

	// Keep additional measurements
	std::vector<double> _xmeasure;
	std::vector<double> _ymeasure;
	int _no_measurements;

};

#endif // !MEASUREMENTDB_H

