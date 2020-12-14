#ifndef SIMULATEDANNEALINGEXTRAP_H
#define SIMULATEDANNEALINGEXTRAP_H

#include "LibraryInterface.h"
#include <map>

SimAnMo::FunctionModel findBestModel(std::map<double, double>& training_points,
	std::map<double, double>& measurement_points, std::string options);

#endif
