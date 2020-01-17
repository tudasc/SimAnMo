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
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <random>
#include <omp.h>
#include <limits>
#include "ExtraPSolution.h"
#include "ExponentialSolution.h"
#include "Configurator.h"

using namespace std;


FILE _iob[] = { *stdin, *stdout, *stderr };

#ifdef _WIN32
extern "C" FILE * __cdecl __iob_func(void)
{
	return _iob;
}
#endif

int no_threads = 1;

template<class SolutionType>
int annealingManager() {
	std::string inputfile = Configurator::getInstance().inputfile;
	SolutionType* sol_per_thread = new SolutionType[no_threads];
	MeasurementDBReader dbreader = MeasurementDBReader();
	MeasurementDB* inputDB = dbreader.readInputFile(inputfile);
	CalcuationInfo<SolutionType> calcinf = CalcuationInfo<SolutionType>();
	unsigned int stepcount = 1;

	double tstart = omp_get_wtime();
	doAnnealing<SolutionType>(inputDB, sol_per_thread, calcinf, stepcount);

	// Prepare the report generation	
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

	calcinf.datapoints = inputDB;
	calcinf.runtime = tduration;
	calcinf.RSScost = abs_min_sol.get_costs();

	std::cout << "Found minimal solution cost: " << abs_min_sol.get_costs()
		<< " in " << stepcount << " steps ("
		<< tduration << " s)"
		<< std::endl;
	abs_min_sol.printModelFunction();

	std::cout << abs_min_sol.printModelFunctionLatex().c_str() << std::endl;
	std::cout << abs_min_sol.printModelFunctionLatexShow().c_str() << std::endl;

	CalcuationInfo<ExtraPSolution> calcinf_log = CalcuationInfo<ExtraPSolution>();
	ExtraPSolution abs_min_sol_log;
	// If it is desired to create a solution based on logarithmic curve fitting
	if (Configurator::getInstance().create_log_exp_model)
	{
		double max_log_range_back = Configurator::getInstance().max_log_range;
		Configurator::getInstance().max_log_range = 0.0;
		stepcount = 1;

		ExtraPSolution* sol_per_thread_log = new ExtraPSolution[no_threads];		
		MeasurementDB* inputDB_log = inputDB->cloneToLog2Version(inputDB);
		//inputDB = inputDB_log;

		doAnnealing<ExtraPSolution>(inputDB_log, sol_per_thread_log, calcinf_log, stepcount, false);
		Configurator::getInstance().max_log_range = max_log_range_back;

		min_cost = std::numeric_limits<double>::max();		
		for (int i = 0; i < no_threads; i++) {
			calcinf_log.sol_per_thread.push_back(sol_per_thread_log[i]);
			if (min_cost > sol_per_thread_log[i].get_costs()) {
				min_cost = sol_per_thread_log[i].get_costs();
				abs_min_sol_log = sol_per_thread_log[i];
				calcinf_log.thread_with_solution = i;
			}
		}
		std::cout << "Found minimal solution cost: " << abs_min_sol_log.get_costs()
			<< " in " << stepcount << " steps ("
			<< tduration << " s)"
			<< std::endl;
		abs_min_sol_log.printModelFunction();
		calcinf.min_sol_log = &abs_min_sol_log;
	}

	// LaTeX Config: Comment out to disable
	LatexPrinter<SolutionType> latprint = LatexPrinter<SolutionType>();
	latprint.printSolution("", &abs_min_sol, inputDB, calcinf);

	delete inputDB;
	return 0;
}


template<class SolutionType>
int doAnnealing(MeasurementDB* inputDB, SolutionType* sol_per_thread, CalcuationInfo<SolutionType>& calcinf,
	unsigned int& stepcount = 1, bool do_quality_log = true) {
	RSSCostCalculator refCostCalc = RSSCostCalculator(inputDB);
#ifdef USE_NAG
	ParameterEstimator paramest = ParameterEstimator(inputDB);
#else
	EigenParameterEstimator paramest = EigenParameterEstimator(inputDB);
#endif

	double ref_array[5] = { 37.760, 0.360, 0.300, 1.0, 2.0 }; // LLLRRDelta

	SolutionType ref_sol = SolutionType(ref_array);
	refCostCalc.calculateCost(&ref_sol);
	cout << "Reference solution cost: " << ref_sol.get_costs() << endl;

	TemperatureInitializer<SolutionType> tempin = TemperatureInitializer<SolutionType>(inputDB);
	double temp_init = tempin.estimateInitialCost(50, 32);
	double T = temp_init;

 	Configurator::getInstance().num_threads = no_threads;

#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		SolutionType act_sol = SolutionType();
		SolutionType abs_min_sol_thread = SolutionType();
		SolutionModifier<SolutionType, RSSCostCalculator> solmod = SolutionModifier<SolutionType, RSSCostCalculator>(inputDB);
		StartSolutionFinder<SolutionType, RSSCostCalculator> startfind = StartSolutionFinder<SolutionType, RSSCostCalculator>(inputDB);

#pragma omp barrier

#pragma omp critical
		startfind.findStartSolution(&act_sol, tid, no_threads);

#pragma omp critical
		{
			std::cout << "Initial cost in thread: " << omp_get_thread_num() << " is "
				<< act_sol.get_costs() << " with temp: " << temp_init << std::endl;
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

		while (T > 0.000001) { // 0.000001
			//cout << "Down";
			for (int i = 0; i < 100; i++) {
				// Generate new solution candidate
				act_sol = solmod.randomModifySolution(&act_sol);

				// Accept solution since it is better
				if (act_sol.get_costs() < min_sol.get_costs()) {
					min_sol = act_sol;
					
					if (do_quality_log) {
						std::pair<unsigned int, double> newpair(stepcount, min_sol.get_costs());
						QualityLogger::getInstance().insertEntry(newpair, tid);
					}

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
						if (do_quality_log) {
							std::pair<unsigned int, double> newpair(stepcount, act_sol.get_costs());
							QualityLogger::getInstance().insertEntry(newpair, tid);
						}
					}
					else
						act_sol = min_sol;
				}
				else
					act_sol = min_sol;
#pragma omp atomic
				stepcount++;
			}
			T = T * 0.9;
		}

		sol_per_thread[omp_get_thread_num()] = abs_min_sol_thread;
	} // End OMP parallel



	// Comparison with linear regression model
#ifdef USE_NAG
	LinearRegressionFinder linfind = LinearRegressionFinder(inputDB);
#else
	AlglibLinearRegressionFinder linfind = AlglibLinearRegressionFinder(inputDB);
#endif

	LinearSolution linsol = linfind.findSolution();
	linsol.printModelFunction();

	calcinf.iterations = stepcount;	
	calcinf.lin_sol = linsol;
	calcinf.print_ref_solution = true;
	calcinf.ref_solution = ref_sol;
	calcinf.print_measurepoints = true;	
	return 0;
}

void printHelp() {
	cout << "Usage: SimulatedAnnealingExtraP [options] "
		<< " --inputfile PATH_AND_NAME_OF_INPUTFILE --outpath PATH_WHERE_TO_PLACE_OUTPUT " << endl
		<< "\t\t" << " --texfile NAME_OF_TEX_AND_PDF_FILES" << endl
		<< endl << endl;
	cout << "--help / -h" << setw(35) << "Print the help" << endl;
	cout << "--number_of_threads / --nt" << setw(35) << "How many threads anneal in parallel (default=1)" << endl;
	cout << "--print_confidence / --pc" << setw(45) << "Print the confidence interval in the predictiion (default=false)" << endl
		<< "--confidence_interval / --ci" << setw(45) << "Set size of confidence interval when printing it (default=0.0)" << endl;
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

		if (input == "--help" || input == "-h") {
			printHelp();
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

		if (input == "--max_log_range" || input == "--logpol") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter max_log_range. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().max_log_range = atof(argv[i + 1]);
			i++;
		}		

		if (input == "--number_of_threads" || input == "--nt") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter number if threads. Terminating." << std::endl;
				exit(-1);
			}
			no_threads = atoi(argv[i + 1]);
			Configurator::getInstance().num_threads = atoi(argv[i + 1]);
			i++;
		}

		// Printing Configuration
		if (input == "--confidence_interval" || input == "--ci") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter epol. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().confidence_interval = atof(argv[i + 1]);
			i++;
		}

		if (input == "--print_confidence" || input == "--pc") {
			Configurator::getInstance().print_confidence = true;
		}

		
	}

	if(Configurator::getInstance().texfile=="")
		Configurator::getInstance().texfile = "DefaultModel";

	if (Configurator::getInstance().outpath == "")
		Configurator::getInstance().outpath = "C:\\temp";



	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(no_threads); // Use X threads for all consecutive parallel regions

	annealingManager<ExponentialSolution>();
	//annealingManager<ExtraPSolution>();
	return 0;

	//MeasurementDB* inputDB = dbreader.giveExampleMeasurementDB05();
	//MeasurementDB* inputDB = dbreader.readInputFile("C:\\temp\\ldsievem3.txt");
	//MeasurementDB* inputDB = dbreader.readInputFile("C:\\temp\\hashsieve.txt");	

	// Comparison to real solution, if available
	/*double refvals[5] = { 0.25, 0.0003,2.25,0.5,0 };
	Solution ref_sol = Solution(refvals);*/

    return 0;
}

