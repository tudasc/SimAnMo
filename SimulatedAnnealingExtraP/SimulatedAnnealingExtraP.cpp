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
#include "GeneralParameterEstimator.h"
#include "AlglibLinearRegressionFinder.h"
#endif


#include "StartSolutionFinder.h"
#include "LatexPrinter.h"
#include "QualityLogger.h"
#include "RingQueue.h"
#include "TemperatureInitializer.h"
#include <iostream>
#include <sstream>      // std::istringstream
#include <iterator>
#include <stdio.h>
#include <math.h>
#include <random>
#include <omp.h>
#include <algorithm>    // std::transform
#include <limits>

#include "ExtraPSolution.h"
#include "ExponentialSolution.h"
#include "ExponentialPolynomSolution.h"
#include "FactorialSolution.h"
#include "Configurator.h"
#include "SimAnMoTools.h"
#include "SimulatedAnnealingExtraP.h"

#include "raRSDParameterEstimator.h"

using namespace std;


FILE _iob[] = { *stdin, *stdout, *stderr };

#ifdef _WIN32
extern "C" FILE * __cdecl __iob_func(void)
{
	return _iob;
}
#endif

//int no_threads = 1;




template<class SolutionType, class CostCalculatorType>
double doAnnealing(MeasurementDB* inputDB, SolutionType* sol_per_thread, CalcuationInfo<SolutionType>& calcinf,
	unsigned int& stepcount, bool do_quality_log = false,
	int steps_per_it = 25, double target_temp = 1e-9, double _cooling_rate = 0.99) {
	CostCalculatorType refCostCalc = CostCalculatorType(inputDB);
#ifdef USE_NAG
	ParameterEstimator paramest = ParameterEstimator(inputDB);
#else
	GeneralParameterEstimator paramest = GeneralParameterEstimator(inputDB);
#endif
	stepcount = 1;

	raRSDParameterEstimator newesti = raRSDParameterEstimator(inputDB);

	//double ref_array[5] = { 25, 3.75E-18, 0.1, 1, 0.0 }; // LLLRRDelta

	//double ref_array[5] = { 2.71575, 3.31285e-09, 1.00153e-06, 0.0, 0.0 }; // BestSolFac
	//double ref_array[5] = { 0.101914 , 0.0713363, -0.146112 , 9.07771e-79 , 0 }; // ManSolFac

	double ref_array[5] = { 0 , 0, 0.745 , 1.9 , 0 };
	SolutionType ref_sol = SolutionType(ref_array);
	//newesti.estimateParameters(&ref_sol);
	//refCostCalc.calculateCost(&ref_sol);
	//exit(123);
	int ret = paramest.estimateParameters(&ref_sol);
	if (ret > 0 && ret < 100) {
		cout << "WARNING: Reference Solution Invalid!";
	}

	else {
		refCostCalc.calculateCost(&ref_sol);
	}


	//cout << "Reference solution cost: " << ref_sol.get_costs() << endl;
	int steps = 0;
	//TemperatureInitializer<SolutionType, CostCalculatorType> tempin = TemperatureInitializer<SolutionType, CostCalculatorType>(inputDB);
	double temp_init = 1.0e-1;// tempin.estimateInitialCost(550, 32);
	//target_temp = temp_init * target_temp;

	const double cooling_rate = _cooling_rate;
	const int step_max = steps_per_it;

	const int no_threads = Configurator::getInstance().num_threads;

#pragma omp parallel num_threads( no_threads )
	{
		int tid = omp_get_thread_num();
		double T =  temp_init;
		SolutionType act_sol = SolutionType();
		SolutionType abs_min_sol_thread = SolutionType();
		SolutionModifier<SolutionType, CostCalculatorType> solmod = SolutionModifier<SolutionType, CostCalculatorType>(inputDB);
		StartSolutionFinder<SolutionType, CostCalculatorType> startfind = StartSolutionFinder<SolutionType, CostCalculatorType>(inputDB);
#pragma omp barrier

//#pragma omp critical
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
		std::uniform_real_distribution<double> distreal(0.0001, 1.0);
#if SIMANMO_VERBOSE > 0
		cout << "Starting temp is: " << target_temp << endl;
#endif

		double progress = 0.0;
#if SIMANMO_VERBOSE > 0
		const int barWidth = 10;
#endif
		const double progressstepwidth = 1.0 / ((log(target_temp/T) / log(cooling_rate)));		

		int without_glob_improve = 0;
		int without_glob_improve2 = 0;		
		while (T > target_temp) {  

#if SIMANMO_VERBOSE > 0
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
#endif
			/*if (act_sol.get_costs() > 10E120Configurator::getInstance().max_cost) {
				cout << act_sol.get_costs() << endl;
				int stop = 1;
				//cin >> stop;
			}*/
			progress += progressstepwidth; // for demonstration only
			steps++;

			for (int i = 0; i < step_max; i++) {

				if (without_glob_improve == 2000) {
					act_sol = abs_min_sol_thread;
					without_glob_improve = 0;
					//cout << "Backtracked";
				}

				if (without_glob_improve2 == 250000) {
					cout << "Thread " << tid << " ends." << endl;
					T = 0;
					break;
				}

				// Generate new solution candidate
				SolutionType act_sol_now = solmod.randomModifySolution(&act_sol);

					if (act_sol_now.get_costs() < 1e-10)
					{
						act_sol_now.printModelFunction();
						cout << "ERRR" << endl;
						int stop = 1;
						//cin >> stop;
					}
				//cout << "c";
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

					//
					double accept_prob = distreal(rng);
					

					if (prob > accept_prob) {
						//cout << prob << " : " << act_sol_now.get_costs() / act_sol.get_costs() << endl;
						if (act_sol_now.get_costs() / act_sol.get_costs() > 2) {
							/*cout << accept_prob << " vs. " << prob 
								<< " with T:" << T << endl;
							int stop = 1;
							//cin >> stop;*/
						}
						//cout << prob << " : " << act_sol_now.get_costs() / act_sol.get_costs() << endl;

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
#if SIMANMO_VERBOSE > 1
	linsol.printModelFunction();
#endif

	refCostCalc.calculateCost(&sol_per_thread[0]);
	cout << "Cost before: " << sol_per_thread[0].get_costs() << endl;
	
	sol_per_thread[0].updateAt(0, 0.101966);
	paramest.estimateParameters(&sol_per_thread[0]);
	cout << "Cost between: " << sol_per_thread[0].get_costs() << endl;

	sol_per_thread[0].updateAt(0, 0.101966);
	refCostCalc.calculateCost(&sol_per_thread[0]);
	cout << "Cost after: " << sol_per_thread[0].get_costs() << endl;

	calcinf.iterations = stepcount;	
	calcinf.lin_sol = linsol;
	calcinf.print_ref_solution = true;
	calcinf.ref_solution = ref_sol;
	calcinf.print_measurepoints = true;	
	return 0;
}

template<class SolutionType, class CostCalcType>
SolutionType annealingManager(MeasurementDB* idb = nullptr) {
	MeasurementDB* inputDB = nullptr;
	if (idb != nullptr) {
		inputDB = idb;
	}

	else {
		if (SimAnMo::readInputMeasurementData(inputDB)) {
			cerr << "Exiting due to measurement data reading error." << endl;
			exit(1001);
		}
	}

	if (Configurator::getInstance().gen_robust) {
		cout << "Will generate robust model." << endl;
		inputDB->unifyMeasurementsToTraining();
	}

	SolutionType* sol_per_thread = new SolutionType[Configurator::getInstance().num_threads ];
	CalcuationInfo<SolutionType> best_calcinf = CalcuationInfo<SolutionType>();
	double best_cost = std::numeric_limits<double>::infinity();
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

		doAnnealing<SolutionType, CostCalcType>(inputDB, sol_per_thread, calcinf, stepcount, Configurator::getInstance().print_costs,
			Configurator::getInstance().ann_steps, Configurator::getInstance().ann_target_temp, Configurator::getInstance().ann_cooling_rate);

		const int no_threads = Configurator::getInstance().num_threads;

		// Prepare the report generation	
		// Get the minimal solution out of all
		min_cost = std::numeric_limits<double>::max();
		SolutionType abs_min_sol;
		for (int i = 0; i < no_threads; i++) {
			calcinf.sol_per_thread.push_back(sol_per_thread[i]);
#if SIMANMO_VERBOSE > 1
			cout << "Found minimal cost in thread " << i << " is " << sol_per_thread[i].get_costs() << " for" << endl;
			sol_per_thread[i].printModelFunction();
#endif

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



		calcinf.runtime = tduration;
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

		ExtraPSolution* sol_per_thread_log = new ExtraPSolution[Configurator::getInstance().num_threads];
		MeasurementDB* inputDB_log = inputDB->cloneToLogVersion(inputDB);
		//inputDB = inputDB_log;

		doAnnealing<ExtraPSolution, RSSCostCalculator>(inputDB_log, sol_per_thread_log, calcinf_log, stepcount, false,
			10, 1e-8, 0.99);

		Configurator::getInstance().max_log_range = max_log_range_back;
		Configurator::getInstance().max_pol_range = max_pol_range_back;
		Configurator::getInstance().min_pol_range = min_pol_range_back;

		min_cost = std::numeric_limits<double>::max();
		for (int i = 0; i < Configurator::getInstance().num_threads; i++) {
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
	
	//delete inputDB;
	return best_abs_min_sol;
}

SimAnMo::FunctionModel findBestModel(std::map<double, double>& training_points,
	std::map<double, double>& measurement_points, std::string options) {

	// Input parameters to standardformat
	istringstream iss(options);
	vector<string> voptions{ istream_iterator<string>{iss},
						  istream_iterator<string>{} };

	cout << "Found options: ";
	for (auto i : voptions) {
		cout << i << " ";
	}
	cout << endl;

	std::vector<char*> cstrings;
	cstrings.reserve(voptions.size());
	for (size_t i = 0; i < voptions.size(); ++i)
		cstrings.push_back(const_cast<char*>(voptions[i].c_str()));

	int temp;
	if (!cstrings.empty())
		SimAnMo::parseConsoleParameters((int)cstrings.size(), &cstrings[0], temp);

	// Input DB from maps
	MeasurementDB mdb = MeasurementDB(training_points, measurement_points);

	// Test a Polynomial Logarithmical Solution
	ExtraPSolution exsol = annealingManager<ExtraPSolution, nnrRSSCostCalculator>(&mdb);
	SimAnMo::FunctionModel funcmod_pollog = SimAnMo::FunctionModel(new ExtraPSolution(exsol));

	// Test a Exponential Solution
	if (funcmod_pollog.getCosts() > 0.1) {
		ExponentialPolynomSolution expolsol =
			annealingManager<ExponentialPolynomSolution, nnrRSSCostCalculator>(&mdb);
		SimAnMo::FunctionModel funcmod_expol =
			SimAnMo::FunctionModel(new ExponentialPolynomSolution(expolsol));

		if (funcmod_pollog.getCosts() < funcmod_expol.getCosts())
			return funcmod_pollog;

		else
			return funcmod_expol;
	}


	return funcmod_pollog;
	
}

int findAModel(std::string mtype, std::string costcaltype) {
	std::transform (mtype.begin(), mtype.end(), mtype.begin(), ::tolower);
	std::transform (costcaltype.begin(), costcaltype.end(), costcaltype.begin(), ::tolower);

	if(mtype.compare("extrapsolution") == 0 && costcaltype.compare("nnrrsscostcalculator")==0) {
		cout << "extrapsolution/nnrrsscostcalculator" << endl;
		//annealingManager<Solution>();
		//annealingManager<ExponentialSolution, nnrRSSCostCalculator>();
		annealingManager<ExtraPSolution, nnrRSSCostCalculator>();
	}

	else if (mtype.compare("extrapsolution") == 0 && costcaltype.compare("rsscostcalculator") == 0) {
		cout << "extrapsolution/nnrrsscostcalculator" << endl;
		//annealingManager<Solution>();
		//annealingManager<ExponentialSolution, nnrRSSCostCalculator>();
		annealingManager<ExtraPSolution, RSSCostCalculator>();
	}

	else if (mtype.compare("factorialsolution") == 0 && costcaltype.compare("nnrrsscostcalculator") == 0) {
		cout << "factorialsolution/nnrrsscostcalculator" << endl;
		annealingManager<FactorialSolution, nnrRSSCostCalculator>();
	}

	else if (mtype.compare("factorialsolution") == 0 && costcaltype.compare("rsscostcalculator") == 0) {
		cout << "factorialsolution/rsscostcalculator" << endl;
		annealingManager<FactorialSolution, RSSCostCalculator>();
	}

	else if (mtype.compare("exponentialsolution") == 0 && costcaltype.compare("nnrrsscostcalculator") == 0) {
		cout << "exponentialsolution/nnrrsscostcalculator" << endl;
		annealingManager<ExponentialPolynomSolution, nnrRSSCostCalculator>();
	}

	else if (mtype.compare("exponentialsolution") == 0 && costcaltype.compare("rsscostcalculator") == 0) {
		cout << "exponentialsolution/rsscostcalculator" << endl;
		annealingManager<ExponentialPolynomSolution, RSSCostCalculator>();
	}

	else {
		cout << "No match for " << mtype << " and " << costcaltype << endl;
	}

	return 0;
}