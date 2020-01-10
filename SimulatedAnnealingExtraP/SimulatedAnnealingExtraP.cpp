#include "MeasurementDB.h"
#include "MeasurementDBReader.h"
#include "Solution.h"
#include "SolutionModifier.h"
#include "RSSCostCalculator.h"
#ifdef USE_NAG
#include "ParameterEstimator.h"
#include "LinearRegressionFinder.h"
#else
#include "EigenParameterEstimator.h"
#include "AlglibLinearRegressionFinder.h"
#endif


#include "StartSolutionFinder.h"
#include "LatexPrinter.h"
#include "QualityLogger.h"
#include "RingQueue.h"
#include "TemperatureInitializer.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <random>
#include <omp.h>
#include <limits>
#include "ExtraPSolution.h"
#include "Configurator.h"


FILE _iob[] = { *stdin, *stdout, *stderr };

#ifdef _WIN32
extern "C" FILE * __cdecl __iob_func(void)
{
	return _iob;
}
#endif

int no_threads = 1;


template<class SolutionType>
int doAnnealing() {
	std::string inputfile = Configurator::getInstance().inputfile;
	SolutionType* sol_per_thread = new SolutionType[no_threads];
	MeasurementDBReader dbreader = MeasurementDBReader();
	MeasurementDB* inputDB = dbreader.readInputFile(inputfile);
	RSSCostCalculator refCostCalc = RSSCostCalculator(inputDB);
#ifdef USE_NAG
	ParameterEstimator paramest = ParameterEstimator(inputDB);
#else
	EigenParameterEstimator paramest = EigenParameterEstimator(inputDB);
#endif
	//double ref_array[4] = { 3.76353, 1.94033e-9, 5, 1.5 };
	//double ref_array[4] = { -30.719, 4.468e-06, 2.5, 2.5 }; // n1.00
	//double ref_array[4] = { 4.454, 2.576e-08, 4.5, 0 }; // n1.50
	//double ref_array[4] = { -7.521, 2.259e-09, 4.5, 2 }; // n1.75
	//double ref_array[4] = { 3.764, 1.940e-09, 5, 1.5 }; // n2.00
	//double ref_array[4] = { -0.962, 4.191e-10, 5, 3 }; // n2.25
	//double ref_array[4] = { -1.336, 3.486e-10, 5, 3.5 }; // n2.40
	//double ref_array[4] = { -0.10, 4.534e-07, 3.5, 1.5 }; // n1.00 NTL
	//double ref_array[4] = { -0.277, 1.526e-07, 3.5, 1.5 }; // n1.00 fplll

	// Exponential
	//double ref_array[5] = { 4.809, 9.594e-10, 1.0, 0.4, 1.0}; // HashSieve full
	//double ref_array[5] = { -0.211, 1.780e-9, 0.350, 1.0, 1.5 }; // HashSieve
	//double ref_array[5] = { -51.973, 6.055e-08, 0.300, 1.0, 1.5 }; // LDSieveM4
	//double ref_array[5] = { -18.080, 4.106e-09, 0.300, 1.0, 2.0 }; // LDSieveM3
	double ref_array[5] = { 37.760, 0.360, 0.300, 1.0, 2.0 }; // LLLRRDelta


	SolutionType ref_sol = SolutionType(ref_array);
	refCostCalc.calculateCost(&ref_sol);
	cout << "Reference solution cost: " << ref_sol.get_costs() << endl;

	TemperatureInitializer<SolutionType> tempin = TemperatureInitializer<SolutionType>(inputDB);
	double temp_init = tempin.estimateInitialCost(1000, 32);
	double T = temp_init;

	double tstart = omp_get_wtime();
	unsigned int stepcount = 1;

	Configurator::getInstance().num_threads = no_threads;

#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		SolutionType act_sol = SolutionType();
		SolutionType abs_min_sol_thread = SolutionType();
		SolutionModifier<SolutionType, RSSCostCalculator> solmod = SolutionModifier<SolutionType, RSSCostCalculator>(inputDB);
		StartSolutionFinder<SolutionType, RSSCostCalculator> startfind = StartSolutionFinder<SolutionType, RSSCostCalculator>(inputDB);

		//cout << tid << " waits." << endl;
#pragma omp barrier
		//cout << tid << " goes on." << endl;

#pragma omp critical
		startfind.findStartSolution(&act_sol, tid, no_threads);

#pragma omp critical
		{
			std::cout << "Initial cost in thread: " << omp_get_thread_num() << " is "
				<< act_sol.get_costs() << std::endl;
			// Step 0 is the initial cost
			std::pair<unsigned int, double> newpair(0, act_sol.get_costs());
			QualityLogger::getInstance().insertEntry(newpair, tid);
		}

		act_sol.printModelFunction();
		SolutionType min_sol = act_sol;
		abs_min_sol_thread = act_sol;

		// For decision if worse solution is accepted
		std::mt19937 rng;
		rng.seed(std::random_device()());
		std::uniform_real_distribution<double> distreal(0.0, 1.0);

		while (T > 000000001) { // 0.000001
			//cout << "Down";
			for (int i = 0; i < 100; i++) {
				// Generate new solution candidate
				act_sol = solmod.randomModifySolution(&act_sol);

				// Accept solution since it is better
				if (act_sol.get_costs() < min_sol.get_costs()) {
					min_sol = act_sol;
					std::pair<unsigned int, double> newpair(stepcount, min_sol.get_costs());
					QualityLogger::getInstance().insertEntry(newpair, tid);

					if (act_sol.get_costs() < abs_min_sol_thread.get_costs()) {
						abs_min_sol_thread = act_sol;
					}
				}

				// Randomize if the worse solution is accepted
				else if (act_sol.get_costs() > min_sol.get_costs()) {
					double prob = exp(-(act_sol.get_costs() - min_sol.get_costs()) / T);
					double accept_prob = distreal(rng);

					if (accept_prob < prob) {
						min_sol = act_sol;
						std::pair<unsigned int, double> newpair(stepcount, act_sol.get_costs());
						QualityLogger::getInstance().insertEntry(newpair, tid);
					}
					else
						act_sol = min_sol;
				}
				else
					act_sol = min_sol;
#pragma omp atomic
				stepcount++;
			}
			T = T * 0.998;
		}

		sol_per_thread[omp_get_thread_num()] = abs_min_sol_thread;
	} // End OMP parallel

	  // Prepare the report generation
	CalcuationInfo<SolutionType> calcinf = CalcuationInfo<SolutionType>();
	// Get the minimal solution out of all
	double min_cost = std::numeric_limits<double>::max();
	SolutionType abs_min_sol;
	for (int i = 0; i < no_threads; i++) {
		calcinf.sol_per_thread.push_back(sol_per_thread[i]);
		if (min_cost > sol_per_thread[i].get_costs()) {
			min_cost = sol_per_thread[i].get_costs();
			abs_min_sol = sol_per_thread[i];
			calcinf.thread_with_solution = i;
		}
	}

	double tduration = omp_get_wtime() - tstart;

	// Comparison with linear regression model
#ifdef USE_NAG
	LinearRegressionFinder linfind = LinearRegressionFinder(inputDB);
#else
	AlglibLinearRegressionFinder linfind = AlglibLinearRegressionFinder(inputDB);
#endif

	LinearSolution linsol = linfind.findSolution();
	linsol.printModelFunction();

	calcinf.iterations = stepcount;
	calcinf.RSScost = abs_min_sol.get_costs();
	calcinf.runtime = tduration;
	calcinf.datapoints = inputDB;
	calcinf.lin_sol = linsol;
	calcinf.print_ref_solution = true;
	calcinf.ref_solution = ref_sol;
	calcinf.print_measurepoints = true;

	std::cout << "Found minimal solution cost: " << abs_min_sol.get_costs()
		<< " in " << stepcount << " steps ("
		<< tduration << " s)"
		<< std::endl;
	abs_min_sol.printModelFunction();

	std::cout << abs_min_sol.printModelFunctionLatex().c_str() << std::endl;
	std::cout << abs_min_sol.printModelFunctionLatexShow().c_str() << std::endl;

	// LaTeX Config: Comment out to disable
	LatexPrinter<SolutionType> latprint = LatexPrinter<SolutionType>();
	latprint.printSolution("", &abs_min_sol, inputDB, calcinf);

	delete inputDB;
	return 0;
}

int main(int argc, char** argv)
{
#ifdef USE_NAG
	cout << "Running Modeler with NAG-Support" << endl;
#endif
	for (int i = 1; i < argc; i++) {
		std::string input = std::string(argv[i]);

		if (input == "--inputfile") {
			if (argc <= i) {
				std::cerr << "Missing inputfile in argument inputfile. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().inputfile = std::string(argv[i + 1]);
			i++;
		}

		if (input == "--texfile") {
			if (argc <= i) {
				std::cerr << "Missing texfile in argument inputfile. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().texfile = std::string(argv[i + 1]);
			i++;
		}

		if (input == "--outpath") {
			if (argc <= i) {
				std::cerr << "Missing inputfile in argument outpath. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().outpath = std::string(argv[i + 1]);
			i++;
		}

		if (input == "--elog") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter elog. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().max_log_range = atof(argv[i + 1]);
			i++;
		}

		if (input == "--epol") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter epol. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().max_pol_range = atof(argv[i + 1]);
			i++;
		}

		if (input == "--nt") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter number if threads. Terminating." << std::endl;
				exit(-1);
			}
			no_threads = atoi(argv[i + 1]);
			Configurator::getInstance().num_threads = atoi(argv[i + 1]);
			i++;
		}

		
	}

	if(Configurator::getInstance().texfile=="")
		Configurator::getInstance().texfile = "DefaultModel";

	if (Configurator::getInstance().outpath == "")
		Configurator::getInstance().outpath = "C:\\temp";



	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(no_threads); // Use X threads for all consecutive parallel regions

	doAnnealing<Solution>();
	return 0;

	//MeasurementDB* inputDB = dbreader.giveExampleMeasurementDB05();
	//MeasurementDB* inputDB = dbreader.readInputFile("C:\\temp\\ldsievem3.txt");
	//MeasurementDB* inputDB = dbreader.readInputFile("C:\\temp\\hashsieve.txt");	

	// Comparison to real solution, if available
	/*double refvals[5] = { 0.25, 0.0003,2.25,0.5,0 };
	Solution ref_sol = Solution(refvals);*/

    return 0;
}

