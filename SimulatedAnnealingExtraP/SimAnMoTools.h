#ifndef  SIMANMOTOOLS_H
#define SIMANMOTOOLS_H


#include "Configurator.h"
#include "MeasurementDB.h"

namespace SimAnMo {

void test();
int parseConsoleParameters(int argc, char** argv, int& no_threads);
void printHelp();
int readInputMeasurementData(MeasurementDB*& inputDB);

}
#endif // ! SIMANMOTOOLS_H

