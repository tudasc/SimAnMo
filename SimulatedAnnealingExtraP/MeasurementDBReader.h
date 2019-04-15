#ifndef MEASUREMENTDBREADER_H
#define MEASUREMENTDBREADER_H

#include "MeasurementDB.h"

#include <string>
using namespace std;


class MeasurementDBReader {
public:
	MeasurementDBReader() { };
	~MeasurementDBReader() {};
	MeasurementDB* readInputFile(string file);


	MeasurementDB* giveExampleMeasurementDB01();
	MeasurementDB* giveExampleMeasurementDB02();
	MeasurementDB* giveExampleMeasurementDB03();
	MeasurementDB* giveExampleMeasurementDB04();
	MeasurementDB* giveExampleMeasurementDB05();
};


#endif
