#include <iostream>
#include <fstream>
#include <omp.h>

#include "../SimulatedAnnealingExtraP/LibraryInterface.h"
#include "../SimulatedAnnealingExtraP/InMaps.h"

int main(int argc, char** argv)
{
	std::ofstream out("out6.txt");
	std::streambuf* coutbuf = std::cout.rdbuf(); //save old buf
	//std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out2.txt!

	//std::map<double, double> mtrai{ {1,3.72}, {2,3.74}, {3,3.72}, {4,3.74}, {5,3.72}, {6,3.74}, {17,3.72}, {45,3.74} };
	std::map<double, double> mmess{  };

	SimAnMo::FunctionModel modi = SimAnMo::findModel(m1, mmess, "--pcd --texfile fplll1.00vTest --outpath ../outputs  --nt 4  --ann_steps_wo_mod 20000 --ann_steps 25 --ann_cooling_rate 0.998 --ann_target_temp 1e-14");

	cout << "I found model " << modi.getModelFunction() << " of type " << modi.getTypeOfModelFunction()
		<< " with RSS " << modi.getRSS() << " and arnRSS " << modi.getraRSD()
		<< ". It is constant: " << modi.isConstant() << endl;

	// Mod 2
	modi = SimAnMo::findModel(m2, mmess, "--texfile fplll1.00vTest --outpath ../outputs  --nt 6  --ann_steps_wo_mod 20000 --ann_steps 30 --ann_cooling_rate 0.998 --ann_target_temp 1e-14");

	cout << "I found model " << modi.getModelFunction() << " of type " << modi.getTypeOfModelFunction()
		<< " with RSS " << modi.getRSS() << " and arnRSS " << modi.getraRSD()
		<< ". It is constant: " << modi.isConstant() << endl;

	// Mod 3
	modi = SimAnMo::findModel(m3, mmess, "--texfile fplll1.00vTest --outpath ../outputs  --nt 6  --ann_steps_wo_mod 20000 --ann_steps 30 --ann_cooling_rate 0.998 --ann_target_temp 1e-14");

	cout << "I found model " << modi.getModelFunction() << " of type " << modi.getTypeOfModelFunction()
		<< " with RSS " << modi.getRSS() << " and arnRSS " << modi.getraRSD()
		<< ". It is constant: " << modi.isConstant() << endl;

	// Mod 4
	modi = SimAnMo::findModel(m4, mmess, "--texfile fplll1.00vTest --outpath ../outputs  --nt 6  --ann_steps_wo_mod 20000 --ann_steps 35 --ann_cooling_rate 0.99 --ann_target_temp 1e-14");

	cout << "I found model " << modi.getModelFunction() << " of type " << modi.getTypeOfModelFunction()
		<< " with RSS " << modi.getRSS() << " and arnRSS " << modi.getraRSD()
		<< ". It is constant: " << modi.isConstant() << endl;

	// Mod 5
	modi = SimAnMo::findModel(m5, mmess, "--logy --texfile fplll1.00vTest --outpath ../outputs  --nt 4  --ann_steps_wo_mod 20000 --ann_steps 25 --ann_cooling_rate 0.998 --ann_target_temp 1e-14");

	cout << "I found model " << modi.getModelFunction() << " of type " << modi.getTypeOfModelFunction()
		<< " with RSS " << modi.getRSS() << " and arnRSS " << modi.getraRSD()
		<< ". It is constant: " << modi.isConstant() << endl;


	std::cout.rdbuf(coutbuf); //reset to standard output again

	return 0;
/*#ifdef USE_NAG
	cout << "Running Modeler with NAG-Support" << endl;
#endif
	int depp = 1;
	SimAnMo::parseConsoleParameters(argc, argv, depp);
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(Configurator::getInstance().num_threads); // Use X threads for all consecutive parallel regions

	//annealingManager<Solution>();
	//annealingManager<ExponentialSolution, nnrRSSCostCalculator>();
	//annealingManager<ExponentialPolynomSolution, nnrRSSCostCalculator>();
	//annealingManager<FactorialSolution, nnrRSSCostCalculator>();
	annealingManager<ExtraPSolution, nnrRSSCostCalculator>();
	return 0;*/
}
