#include <iostream>
#include <fstream>
#include <omp.h>

#include "../SimulatedAnnealingExtraP/SimulatedAnnealingExtraP.h"
#include "../SimulatedAnnealingExtraP/SimAnMoTools.h"

#include "../SimulatedAnnealingExtraP/ExtraPSolution.h"
#include "../SimulatedAnnealingExtraP/FactorialSolution.h"
#include "../SimulatedAnnealingExtraP/ExponentialPolynomSolution.h"

#include "../SimulatedAnnealingExtraP/alglib/optimization.h"

using namespace alglib;

int main(int argc, char** argv)
{
#ifdef USE_NAG
	cout << "Running Modeler with NAG-Support" << endl;
#endif
	int depp = 1;
	Configurator::getInstance().param_est_typ = TYPE_ALGLIBPARAMETER;
	SimAnMo::parseConsoleParameters(argc, argv, depp);

	if (Configurator::getInstance().costcalc_type == "rsscostcalculator") {
		cout << "NAM" << endl;
		Configurator::getInstance().param_est_typ = TYPE_EIGENPARAMETER;
	}

	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(Configurator::getInstance().num_threads); // Use X threads for all consecutive parallel regions

	//SCNTContruction
	//double coeff[5] = { 4.3089 , 0.275642 , 0.892573 , 0.992706 , 0 };

	//SCNTContructionPi
	//double coeff[5] = { -0.148519 , 0.365812 , 1.01371 , 0.0019119 , 0 };

	//MMM
	//double coeff[5] = { 0.212195 , 2.86256e-09 , 3.13639 , 0.00345831 , 0 };

	//MMMPi
	//double coeff[5] = { 0.126544 , 9.81281e-07 , 2.9874 , 0.00923343 , 0 };

	//String Permutation
	//double coeff[5] = { 0.101966 , 0.0713363 , -0.146112 , 5.6157e-55 , 0 };

	//String PermutationPi
	//double coeff[5] = { 7.87598 , 0.0423795 , -1.99994 , 3.85656e-42 , 0 };

	//TSP Serial
	//double coeff[5] = { -4.80663e-06 , 1.77276e-08 , -0.432156 , 0 };

	//Test
	double coeff[5] = { 0.101914 , 0.0713363, -0.146112 , 4.32428e-66 , 0 };

	//FactorialSolution exsol = calculateFullModelMetrics<FactorialSolution>(coeff);
	//cout << "RSS: " << exsol._RSS << endl;
	//cout << "nnrRSS: " << exsol._nnrRSS << endl;

	//return 0;
	//return findAModel("extrapsolution", "nnrrsscostcalculator");
	
	//return findAModel("factorialsolution", "nnrrsscostcalculator");
	//return findAModel("factorialsolution", "rsscostcalculator");

	//return findAModel("exponentialsolution", "nnrrsscostcalculator");
	return findAModel(Configurator::getInstance().solution_type, Configurator::getInstance().costcalc_type);
}
