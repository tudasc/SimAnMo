#include <iostream>
#include <fstream>
#include <omp.h>

#include "../SimulatedAnnealingExtraP/SimulatedAnnealingExtraP.h"
#include "../SimulatedAnnealingExtraP/SimAnMoTools.h"



int main(int argc, char** argv)
{
#ifdef USE_NAG
	cout << "Running Modeler with NAG-Support" << endl;
#endif
	int depp = 1;
	SimAnMo::parseConsoleParameters(argc, argv, depp);
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(Configurator::getInstance().num_threads); // Use X threads for all consecutive parallel regions

	std::cout << "I am here " << std::endl;
	//return findAModel("extrapsolution", "nnrrsscostcalculator");
	return findAModel("factorialsolution", "nnrrsscostcalculator");

	//return findAModel("exponentialsolution", "nnrrsscostcalculator");
}
