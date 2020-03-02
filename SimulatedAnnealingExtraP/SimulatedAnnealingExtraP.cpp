#include "MeasurementDB.h"
#include "MeasurementDBReader.h"
#include "Solution.h"
#include "SolutionModifier.h"
#include "RSSCostCalculator.h"
#include "nnrRSSCostCalculator.h"
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
#include "ExponentialPolynomSolution.h"
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


template<class SolutionType, class CostCalculatorType>
double doAnnealing(MeasurementDB* inputDB, SolutionType* sol_per_thread, CalcuationInfo<SolutionType>& calcinf, double target_temp,
	unsigned int& stepcount = 1, bool do_quality_log = false) {
	CostCalculatorType refCostCalc = CostCalculatorType(inputDB);
#ifdef USE_NAG
	ParameterEstimator paramest = ParameterEstimator(inputDB);
#else
	EigenParameterEstimator paramest = EigenParameterEstimator(inputDB);
#endif

	double ref_array[5] = { 25, 3.75E-18, 0.1, 1, 0.0 }; // LLLRRDelta

	SolutionType ref_sol = SolutionType(ref_array);
	refCostCalc.calculateCost(&ref_sol);
	//cout << "Reference solution cost: " << ref_sol.get_costs() << endl;

	TemperatureInitializer<SolutionType, CostCalculatorType> tempin = TemperatureInitializer<SolutionType, CostCalculatorType>(inputDB);
	double temp_init = tempin.estimateInitialCost(250, 32);	

 	Configurator::getInstance().num_threads = no_threads;

#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		double T = temp_init;
		SolutionType act_sol = SolutionType();
		SolutionType abs_min_sol_thread = SolutionType();
		SolutionModifier<SolutionType, CostCalculatorType> solmod = SolutionModifier<SolutionType, CostCalculatorType>(inputDB);
		StartSolutionFinder<SolutionType, CostCalculatorType> startfind = StartSolutionFinder<SolutionType, CostCalculatorType>(inputDB);

#pragma omp barrier

#pragma omp critical
		startfind.findStartSolution(&act_sol, tid, no_threads);

#pragma omp critical
		{
			//std::cout << "Initial cost in thread: " << omp_get_thread_num() << " is "
			//	<< act_sol.get_costs() << " with temp: " << temp_init << std::endl;
			// Step 0 is the initial cost
			std::pair<unsigned int, double> newpair(0, act_sol.get_costs());
			if(do_quality_log)
				QualityLogger::getInstance().insertEntry(newpair, tid);
		}

		//act_sol.printModelFunction();
		//SolutionType min_sol = act_sol;
		abs_min_sol_thread = act_sol;

		// For decision if worse solution is accepted
		std::mt19937 rng;
		rng.seed(std::random_device()());
		std::uniform_real_distribution<double> distreal(0.0, 1.0);

		int without_glob_improve = 0;
		int without_glob_improve2 = 0;

		target_temp = T * target_temp;

		while (T > target_temp) {
			for (int i = 0; i < 150; i++) {

				/*if (without_glob_improve == 15000) {
					act_sol = abs_min_sol_thread;
					without_glob_improve = 0;
				}*/

				if (without_glob_improve2 == 50000) {
					//cout << "Thread " << tid << " ends." << endl;
					T = 0;
					break;
				}

				// Generate new solution candidate
				SolutionType act_sol_now = solmod.randomModifySolution(&act_sol);

				// Accept solution since it is better
				if (act_sol_now.get_costs() < act_sol.get_costs()) {
					act_sol = act_sol_now;

					if (do_quality_log) {
						std::pair<unsigned int, double> newpair(stepcount, act_sol.get_costs());
						QualityLogger::getInstance().insertEntry(newpair, tid);
					}

					if (act_sol.get_costs() < abs_min_sol_thread.get_costs()) {
						abs_min_sol_thread = act_sol;
					}
					else {
						without_glob_improve++;
						without_glob_improve2++;
					}
				}

				// Randomize if the worse solution is accepted
				else if (act_sol_now.get_costs() >= act_sol.get_costs()) {
					double prob = exp(-(act_sol_now.get_costs() - act_sol.get_costs()) / T);
					double accept_prob = distreal(rng);

					//cout << accept_prob  << " :" <<  act_sol_now.get_costs() << " " << act_sol.get_costs() << endl;

					if (accept_prob > prob) {
						act_sol = act_sol_now;
						if (do_quality_log) {
							std::pair<unsigned int, double> newpair(stepcount, act_sol.get_costs());
							QualityLogger::getInstance().insertEntry(newpair, tid);							
						}
					}
					/*else
						act_sol = min_sol;*/
				}
				else {
					cerr << "Illegal flow in Annealing" << endl;   //act_sol = min_sol;
					cerr << act_sol_now.get_costs() << " vs " << act_sol.get_costs() << endl;
					act_sol_now.printModelFunction();
					act_sol.printModelFunction();
					cerr << endl;
					exit(200);
				}
					
#pragma omp atomic
				stepcount++;
			}
			T = T * 0.999;
			//cout << "Target: "<< target_temp << " / T: " << T << endl;
		}

		sol_per_thread[tid] = abs_min_sol_thread;
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

template<class SolutionType>
int annealingManager() {

	Configurator::getInstance().noLogModel();

	std::string inputfile = Configurator::getInstance().inputfile;
	SolutionType* sol_per_thread = new SolutionType[no_threads];
	MeasurementDBReader dbreader = MeasurementDBReader();
	MeasurementDB* inputDB = dbreader.readInputFile(inputfile);
	CalcuationInfo<SolutionType> best_calcinf = CalcuationInfo<SolutionType>();
	double best_cost = std::numeric_limits<double>::max();
	SolutionType best_abs_min_sol;
	double tstart = omp_get_wtime();

	unsigned int stepcount = 1;
	double min_cost = std::numeric_limits<double>::max();
	for (int i = 0; i < 1; i++)
	{
		CalcuationInfo<SolutionType> calcinf = CalcuationInfo<SolutionType>();
		stepcount = 1;
		doAnnealing<SolutionType, nnrRSSCostCalculator>(inputDB, sol_per_thread, calcinf, 1e-11, stepcount, true);

		// Prepare the report generation	
		// Get the minimal solution out of all
		min_cost = std::numeric_limits<double>::max();
		SolutionType abs_min_sol;
		for (int i = 0; i < no_threads; i++) {
			calcinf.sol_per_thread.push_back(sol_per_thread[i]);
			if (min_cost > sol_per_thread[i].get_costs()) {
				min_cost = sol_per_thread[i].get_costs();
				abs_min_sol = sol_per_thread[i];
				calcinf.thread_with_solution = i;
			}
		}

		calcinf.RSScost = abs_min_sol.get_costs();
		double tduration = omp_get_wtime() - tstart;
		std::cout << "Found minimal solution cost: " << abs_min_sol.get_costs()
			<< " in " << stepcount << " steps ("
			<< tduration << " s)"
			<< std::endl;

		best_calcinf.runtime = tduration;

		abs_min_sol.printModelFunction();

		if (min_cost < best_cost) {
			best_cost = min_cost;
			best_calcinf = calcinf;
			best_abs_min_sol = abs_min_sol;
			cout << "Reducing global costs to " << min_cost << endl;
		}
	}	

	best_calcinf.datapoints = inputDB;
	

	CalcuationInfo<ExtraPSolution> calcinf_log = CalcuationInfo<ExtraPSolution>();
	ExtraPSolution abs_min_sol_log;
	// If it is desired to create a solution based on logarithmic curve fitting
	if (Configurator::getInstance().create_log_exp_model)
	{
		double max_log_range_back = Configurator::getInstance().max_log_range;
		double max_pol_range_back = Configurator::getInstance().max_pol_range;
		double min_pol_range_back = Configurator::getInstance().min_pol_range;

		Configurator::getInstance().max_log_range = 0.0;
		Configurator::getInstance().max_pol_range = 1.001;
		Configurator::getInstance().min_pol_range = 0.999;

		stepcount = 1;

		ExtraPSolution* sol_per_thread_log = new ExtraPSolution[no_threads];
		MeasurementDB* inputDB_log = inputDB->cloneToLog2Version(inputDB);
		//inputDB = inputDB_log;

		doAnnealing<ExtraPSolution, RSSCostCalculator>(inputDB_log, sol_per_thread_log, calcinf_log, 0.0005, stepcount, false);
		Configurator::getInstance().max_log_range = max_log_range_back;
		Configurator::getInstance().max_pol_range = max_pol_range_back;
		Configurator::getInstance().min_pol_range = min_pol_range_back;

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
			//<< tduration << " s)"
			<< std::endl;
		abs_min_sol_log.printModelFunction();

		// Replace costs by costs for the non-logarithmized data
		abs_min_sol_log.switchtoLinScale();
		RSSCostCalculator refCostCalc = RSSCostCalculator(inputDB);
		refCostCalc.calculateCost(&abs_min_sol_log);

		best_calcinf.min_sol_log = &abs_min_sol_log;
	}

	// LaTeX Config: Comment out to disable
	LatexPrinter<SolutionType> latprint = LatexPrinter<SolutionType>();
	latprint.printSolution("", &best_abs_min_sol, inputDB, best_calcinf);

	delete inputDB;
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

	//annealingManager<Solution>();
	//annealingManager<ExponentialSolution>();
	annealingManager<ExponentialPolynomSolution>();
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

