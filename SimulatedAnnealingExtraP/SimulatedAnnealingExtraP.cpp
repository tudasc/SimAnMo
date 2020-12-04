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
#include "FactorialSolution.h"
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
double doAnnealing(MeasurementDB* inputDB, SolutionType* sol_per_thread, CalcuationInfo<SolutionType>& calcinf,
	unsigned int& stepcount, bool do_quality_log = false,
	int steps_per_it = 25, double target_temp = 1e-9, double _cooling_rate = 0.99) {
	CostCalculatorType refCostCalc = CostCalculatorType(inputDB);
#ifdef USE_NAG
	ParameterEstimator paramest = ParameterEstimator(inputDB);
#else
	EigenParameterEstimator paramest = EigenParameterEstimator(inputDB);
#endif
	stepcount = 1;
	//double ref_array[5] = { 25, 3.75E-18, 0.1, 1, 0.0 }; // LLLRRDelta

	//double ref_array[5] = { 2.71575, 3.31285e-09, 1.00153e-06, 0.0, 0.0 }; // BestSolFac
	double ref_array[5] = { 0.0, 1, 0.0, 0.0, 0.0 }; // ManSolFac

	SolutionType ref_sol = SolutionType(ref_array);
	refCostCalc.calculateCost(&ref_sol);
	//cout << "Reference solution cost: " << ref_sol.get_costs() << endl;
	int steps = 0;
	//TemperatureInitializer<SolutionType, CostCalculatorType> tempin = TemperatureInitializer<SolutionType, CostCalculatorType>(inputDB);
	double temp_init = 1.0;// tempin.estimateInitialCost(550, 32);
	//target_temp = temp_init * target_temp;

 	Configurator::getInstance().num_threads = no_threads;

    

	const double cooling_rate = _cooling_rate;
	const int step_max = steps_per_it;

#pragma omp parallel num_threads( no_threads )
	{
		int tid = omp_get_thread_num();
		double T =  temp_init;
		SolutionType act_sol = SolutionType();
		SolutionType abs_min_sol_thread = SolutionType();
		SolutionModifier<SolutionType, CostCalculatorType> solmod = SolutionModifier<SolutionType, CostCalculatorType>(inputDB);
		StartSolutionFinder<SolutionType, CostCalculatorType> startfind = StartSolutionFinder<SolutionType, CostCalculatorType>(inputDB);

#pragma omp barrier

#pragma omp critical
		startfind.findStartSolution(&act_sol, tid, no_threads);
		//act_sol = ref_sol;

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
		cout << "Starting temp is: " << target_temp << endl;

		double progress = 0.0;
		const int barWidth = 10;
		const double progressstepwidth = 1.0 / ((log(target_temp/T) / log(cooling_rate)));		

		int without_glob_improve = 0;
		int without_glob_improve2 = 0;		
		while (T > target_temp) {    

			if (tid == 0) {
				std::cout << "[";
				//int pos = barWidth * progress;
				for (int i = 0; i < barWidth; ++i) {
					//if (i < pos) std::cout << "=";
					//else if (i == pos) std::cout << ">";
					//else std::cout << " ";
					//std::cout << " ";
				}
				std::cout << "] " << int(progress * 100.0) << " %\r";
				std::cout.flush();
			}

			progress += progressstepwidth; // for demonstration only
			steps++;

			for (int i = 0; i < step_max; i++) {

				/*if (without_glob_improve == 150000) {
					act_sol = abs_min_sol_thread;
					without_glob_improve = 0;
					cout << "Backtracked";
				}*/

				if (without_glob_improve2 == 1000000) {
					cout << "Thread " << tid << " ends." << endl;
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

					if (accept_prob > prob && 1 == 2) {
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
			T = T * cooling_rate;
			//cout << "Target: "<< target_temp << " / T: " << T << endl;
		}

		sol_per_thread[tid] = abs_min_sol_thread;
	} // End OMP parallel

	cout << steps << " steps." << endl;

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

template<class SolutionType, class CostCalcType>
int annealingManager() {
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
	for (int i = 0; i < Configurator::getInstance().no_of_trials; i++)
	{
		CalcuationInfo<SolutionType> calcinf = CalcuationInfo<SolutionType>();
		stepcount = 1;

		cout << "Annealing with " << Configurator::getInstance().ann_steps << " / " << Configurator::getInstance().ann_target_temp <<
			" / " << Configurator::getInstance().ann_cooling_rate << endl;

		doAnnealing<SolutionType, CostCalcType>(inputDB, sol_per_thread, calcinf, stepcount, true,
			Configurator::getInstance().ann_steps, Configurator::getInstance().ann_target_temp, Configurator::getInstance().ann_cooling_rate);

		// Prepare the report generation	
		// Get the minimal solution out of all
		min_cost = std::numeric_limits<double>::max();
		SolutionType abs_min_sol;
		for (int i = 0; i < no_threads; i++) {
			calcinf.sol_per_thread.push_back(sol_per_thread[i]);

			cout << "Found minimal cost in thread " << i << " is " << sol_per_thread[i].get_costs() << " for" << endl;
			sol_per_thread[i].printModelFunction();

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
		MeasurementDB* inputDB_log = inputDB->cloneToLogVersion(inputDB);
		//inputDB = inputDB_log;

		doAnnealing<ExtraPSolution, RSSCostCalculator>(inputDB_log, sol_per_thread_log, calcinf_log, stepcount, false,
			10, 1e-8, 0.99);

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
		std::cout << "Found minimal solution cost for linear: " << abs_min_sol_log.get_costs()
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
	if (Configurator::getInstance().do_latex_output) {
		LatexPrinter<SolutionType> latprint = LatexPrinter<SolutionType>();
		latprint.printSolution("", &best_abs_min_sol, inputDB, best_calcinf);
	}


	delete inputDB;
	return 0;
}

void printHelp() {
	cout << "Usage: SimulatedAnnealingExtraP [options] "
		<< " --inputfile PATH_AND_NAME_OF_INPUTFILE --outpath PATH_WHERE_TO_PLACE_OUTPUT " << endl
		<< "\t\t" << " --texfile NAME_OF_TEX_AND_PDF_FILES" << endl
		<< endl << endl;
	cout << "--help / -h"  << setw(55) << "Print the help" << endl;
	cout << "" << endl;
	cout << "--inputfile / -i + PATH_TO_INPUT_FILE/FILNAME" << setw(55) << "MANDATORY. The path and name to/of the input file" << endl;
	cout << "--outfile / -o + PATH_TO_OUT_FILE" << setw(55) << "The path to the output files (default=folder of executable)" << endl;
	cout << "--texfile / -t + NAME_OF_FILES" << setw(55) << "MANDATORY. The name of output tex and pdf files" << endl;

	cout << "--number_of_threads / --nt + INT" << setw(55) << "How many threads anneal in parallel (default=1)" << endl;
	cout << "--number_of_trials / --tr + INT" << setw(55) << "How many repetitions of annealing process (default=1)" << endl;	

	cout << "--ann_target_temp / --att + FLOAT" << setw(55) << "Target temperature for annealing (default=1e-9)" << endl;
	cout << "--ann_cooling_rate / --acr + FLOAT" << setw(55) << "Temperature degredation per iteration (default=0.99)" << endl;
	cout << "--ann_steps / --as + INT" << setw(55) << "How many steps are performed per temperature (default=15)" << endl;
	
	cout << "--ann_steps_wo_mod / --awm + INT" << setw(55) << "Heuristic: Stop if no improvement after this number of steps (default=200000)" << endl;
	cout << "--ann_steps_backtrack / --abt + INT" << setw(55) << "Heuristic: Backtrack per thread after this number of steps (default=20000)" << endl;
	

	cout << "--confidence_interval / --ci + FLOAT" << setfill(' ') << setw(55) << "Set size of confidence interval when printing it (default=0.0)" << endl;
	
	// Printing 
	cout << endl << "SECTION: LaTeX Configuration" << endl;
	cout << "--genlatex / --gl" << setfill(' ') << setw(55) << "Activate the LaTeX report generation (default=false)" << endl;
	cout << "--openpdf / --op" << setfill(' ') << setw(55) << "Generated pdf file is automatically opened with pdfxchange at --pathtopdfxchange (default=false)" << endl;
	cout << "--pathtopdfxchange" << setfill(' ') << setw(55) << "If --openpdf is set, pdf file is automatically opened with pdfxchange at this path" << endl;
	cout << "--logy" << setfill(' ') << setw(55) << "The y-axis in the prediction graph is scaled logarithmically" << endl;
	cout << "--print_confidence / --pc" << setw(55) << "Print the confidence interval in the predictiion (default=false)" << endl;
	cout << "--print_cost_details / --pcd" << setw(55) << "Print details of cost development during annealing (default=false)" << endl;

	// Pol-Log
	cout << endl << "SECTION: Polynomial-logarithmic (pol-log) model configuration" << endl;
	cout << "--max_log_range / --melog + FLOAT" << setfill(' ') << setw(55) << "Maximum exponent for lorarithms (default=4.00)" << endl;
	cout << "--min_pol_range / --mipol + FLOAT" << setfill(' ') << setw(55) << "Minimum exponent for polynoms (default=-2.00)" << endl;
	cout << "--max_pol_range / --mepol + FLOAT" << setfill(' ') << setw(55) << "Maximum exponent for polynoms (default=6.00)" << endl;


	// Lin-Log
	cout << endl << "SECTION: Linear-logarithmic (lin-log) model configuration" << endl;
	cout << "--create_lin_log / --ll"  << setfill(' ') << setw(55) << "Create a lin-log model if set (default=false)" << endl;
	cout << "--base_lin_log / --bll + INT" << setfill(' ') << setw(55) << "Basis that is used for lin-log-model (default=2)" << endl;

	// Lin-Log
	cout << endl << "SECTION: Extended Extra-P model configuration" << endl;
	cout << "--max_exp_range / --meexp + INT" << setfill(' ') << setw(55) << "Maximum coefficient in exponent (default=4.00)" << endl;

	// Exp-Pol
	cout << endl << "SECTION: Exponential-polynomial (exp-pol) model configuration" << endl;
	cout << "--exp_pol_min_coeff / --epmic + FLOAT" << setfill(' ') << setw(55) << "Minimum coefficient in the exponent of exp-pol models (default=0.01)" << endl;
	cout << "--exp_pol_max_coeff / --epmac + FLOAT" << setfill(' ') << setw(55) << "Maximum coefficient in the exponent of exp-pol models (default=2.0)" << endl;
	cout << "--exp_pol_min_exp / --epmie + FLOAT" << setfill(' ') << setw(55) << "Minimum exponent in the exponent of exp-pol models (default=0.5)" << endl;
	cout << "--exp_pol_max_exp / --epmae + FLOAT" << setfill(' ') << setw(55) << "Maximum exponent in the exponent of exp-pol models (default=3.0)" << endl;
	exit(0);
}

int main(int argc, char** argv)
{

#ifdef USE_NAG
	cout << "Running Modeler with NAG-Support" << endl;
#endif
	for (int i = 1; i < argc; i++) {
		std::string input = std::string(argv[i]);

		if (input == "--number_of_trials" || input == "--tr") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter number of trials for annealing. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().no_of_trials = atoi(argv[i + 1]);
			i++;
		}

		if (input == "--ann_steps" || input == "--as") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter number of steps int annealing. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().ann_steps = atoi(argv[i + 1]);
			i++;
		}

		if (input == "--ann_cooling_rate" || input == "--acr") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter number of cooling rate for annealing. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().ann_cooling_rate = atof(argv[i + 1]);
			i++;
		}

		if (input == "--ann_target_temp" || input == "--att") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter target temperature for annealing. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().ann_target_temp = atof(argv[i + 1]);
			i++;
		}

		if (input == "--ann_steps_wo_mod" || input == "--awm") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter number of steps without modification for annealing. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().ann_steps_wo_mod = atoi(argv[i + 1]);
			i++;
		}

		if (input == "--ann_steps_backtrack" || input == "--abt") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter number of steps before backtrack for annealing. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().ann_steps_backtrack = atoi(argv[i + 1]);
			i++;
		}
		
		if (input == "-i" || input == "--inputfile") {
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
		

		if (input == "--number_of_threads" || input == "--nt") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter number of threads. Terminating." << std::endl;
				exit(-1);
			}
			no_threads = atoi(argv[i + 1]);
			cout << "Anneal with " << no_threads << " threads." << endl;
			Configurator::getInstance().num_threads = atoi(argv[i + 1]);
			i++;
		}

		if (input == "--number_of_trials" || input == "--tr") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter number of trials for annealing. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().no_of_trials = atoi(argv[i + 1]);
			i++;
		}

		// Printing Configuration
		if (input == "--genlatex" || input == "--gl") {
			Configurator::getInstance().do_latex_output = true;
		}

		if (input == "--openpdf" || input == "--op") {
			Configurator::getInstance().open_latex_output = true;
		}

		if (input == "--pathtopdfxchange") {
			if (argc <= i) {
				std::cerr << "Missing path in argument pathtopdfxchange. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().path_pdf_xchange = std::string(argv[i + 1]);
			i++;
		}

		if (input == "--logy" || input == "--logy") {
			Configurator::getInstance().ymode_log = true;
			cout << "Printing logged" << endl;
		}

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

		if (input == "--print_cost_details" || input == "--pcd") {
			Configurator::getInstance().print_costs = true;
		}


		// Pol-Log model
		if (input == "--max_log_range" || input == "--melog") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter max_log_range. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().max_log_range = atof(argv[i + 1]);
			i++;
		}

		if (input == "--max_pol_range" || input == "--mepol") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter max_pol_range. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().max_pol_range = atof(argv[i + 1]);
			i++;
		}

		if (input == "--min_pol_range" || input == "--mipol") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter mix_pol_range. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().min_pol_range = atof(argv[i + 1]);
			i++;
		}

		// Extended model
		if (input == "--max_exp_range" || input == "--meexp") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter max_exp_range. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().max_exp_range = atof(argv[i + 1]);
			i++;
		}

		// Exp-Pol model
		if (input == "--exp_pol_min_coeff" || input == "--epmic") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter exp_pol-min_coeff. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().min_exp_coeff_range = atof(argv[i + 1]);
			i++;
		}

		if (input == "--exp_pol_max_coeff" || input == "--epmac") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter exp_pol-max_coeff. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().max_exp_coeff_range = atof(argv[i + 1]);
			i++;
		}

		if (input == "--exp_pol_min_exp" || input == "--epmie") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter exp_pol-min_exp. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().min_exp_exp_range = atof(argv[i + 1]);
			i++;
		}

		if (input == "--exp_pol_max_exp" || input == "--epmae") {
			if (argc <= i) {
				std::cerr << "Missing argument for parameter exp_pol-max_exp. Terminating." << std::endl;
				exit(-1);
			}
			Configurator::getInstance().max_exp_exp_range = atof(argv[i + 1]);
			i++;
		}

		// Lin-log Model
		if (input == "--create_lin_log" || input == "--ll") {
			Configurator::getInstance().create_log_exp_model = true;
		}

		if (input == "--base_lin_log" || input == "--bll") {
			if (argc <= i) {
				std::cerr << "Missing argument for base of lin-log model. Terminating." << std::endl;
				exit(-1);
			}			
			Configurator::getInstance().base_for_lin_log = atoi(argv[i + 1]);
			i++;
		}
	}

	if(Configurator::getInstance().texfile=="")
		Configurator::getInstance().texfile = "DefaultModel";

	if (Configurator::getInstance().outpath == "")
		Configurator::getInstance().outpath = "C:\\temp";



	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(no_threads); // Use X threads for all consecutive parallel regions

	//annealingManager<Solution>();
	//annealingManager<ExponentialSolution, nnrRSSCostCalculator>();
	//annealingManager<ExponentialPolynomSolution, nnrRSSCostCalculator>();
	//annealingManager<FactorialSolution, nnrRSSCostCalculator>();
	annealingManager<ExtraPSolution,nnrRSSCostCalculator>();
	return 0;
}

