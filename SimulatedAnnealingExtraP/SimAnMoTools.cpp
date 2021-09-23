#include "SimAnMoTools.h"
#include <omp.h>
#include <iomanip>
#include <iostream>

#include "MeasurementDBReader.h"

using namespace std;

//namespace SimAnMo {

	void SimAnMo::test() {
		cout << "Do something." << endl;
	}

	void SimAnMo::printHelp() {
		cout << "Usage: SimulatedAnnealingExtraP [options] "
			<< " --inputfile PATH_AND_NAME_OF_INPUTFILE --outpath PATH_WHERE_TO_PLACE_OUTPUT " << endl
			<< "\t\t" << " --texfile NAME_OF_TEX_AND_PDF_FILES" << endl
			<< endl << endl;
		cout << "--help / -h" << setw(55) << "Print the help" << endl;
		cout << "" << endl;
		cout << "--inputfile / -i + PATH_TO_INPUT_FILE/FILNAME" << setw(55) << "MANDATORY. The path and name to/of the input file" << endl;
		cout << "--outpath / -o + PATH_TO_OUT_FILE" << setw(55) << "The path to the output files (default=folder of executable)" << endl;
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
		cout << "--min_log_range / --milog + FLOAT" << setfill(' ') << setw(55) << "Minimum exponent for lorarithms (default=0.00)" << endl;
		cout << "--max_log_range / --melog + FLOAT" << setfill(' ') << setw(55) << "Maximum exponent for lorarithms (default=4.00)" << endl;
		cout << "--min_pol_range / --mipol + FLOAT" << setfill(' ') << setw(55) << "Minimum exponent for polynoms (default=-2.00)" << endl;
		cout << "--max_pol_range / --mepol + FLOAT" << setfill(' ') << setw(55) << "Maximum exponent for polynoms (default=6.00)" << endl;


		// Lin-Log
		cout << endl << "SECTION: Linear-logarithmic (lin-log) model configuration" << endl;
		cout << "--create_lin_log / --ll" << setfill(' ') << setw(55) << "Create a lin-log model if set (default=false)" << endl;
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
		
		// Fac-Pol
		cout << endl << "SECTION: Factorial (fac) model configuration" << endl;
		cout << "--fac_pol_min_coeff / --fpmic + FLOAT" << setfill(' ') << setw(55) << "Minimum coefficient in the exponent of the polynomial part of the fac models (default=-0.5)" << endl;
		cout << "--fac_pol_max_coeff / --fpmac + FLOAT" << setfill(' ') << setw(55) << "Maximum coefficient in the exponent of the polynomial part of the fac models (default=3.0)" << endl;
		cout << "--fac_log_min_coeff / --flmic + FLOAT" << setfill(' ') << setw(55) << "Minimum coefficient in the exponent of the logarithmic part of the fac models (default=1e-3)" << endl;
		cout << "--fac_log_max_coeff / --flmac + FLOAT" << setfill(' ') << setw(55) << "Maximum coefficient in the exponent of the logarithmic part of the fac models (default=1.5)" << endl;

		exit(0);
	}

	int SimAnMo::parseConsoleParameters(int argc, char** argv, int& no_threads) {
		for (int i = 0; i < argc; i++) {
			std::string input = std::string(argv[i]);
			cout << "Parsing: " << input << endl;

			if (input == "--solution_type" || input == "--st") {
				if (argc <= i) {
					std::cerr << "Missing argument for parameter solution type. Terminating." << std::endl;
					exit(-1);
				}
				Configurator::getInstance().solution_type = string(argv[i + 1]);
				i++;
			}

			if (input == "--costcalc_type" || input == "--cct") {
				if (argc <= i) {
					std::cerr << "Missing argument for parameter cost calcuator type. Terminating." << std::endl;
					exit(-1);
				}
				Configurator::getInstance().costcalc_type = string(argv[i + 1]);
				i++;
			}

			if (input == "--paramest_type" || input == "--pet") {
				if (argc <= i) {
					std::cerr << "Missing argument for parameter estimation type. Terminating." << std::endl;
					exit(-1);
				}
				string strtoparse = string(argv[i + 1]);

				if (strtoparse.compare("rarsdest") == 0 &&
					(Configurator::getInstance().costcalc_type.compare("rarsdcost")==0
					|| Configurator::getInstance().costcalc_type.compare("nnrrsscostcalculator") == 0
						)) {
					Configurator::getInstance().param_est_typ = TYPE_EIGENPARAMETER;
					cout << "Choosing experimentally raRSD parameter estimator"
						<< endl;
				}
					
				else
					Configurator::getInstance().param_est_typ = TYPE_EIGENPARAMETER;

				Configurator::getInstance().param_est_typ = 
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
				std::cout << "Outpath is: " << input << std::endl;
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

			if (input == "--generate_robust" || input == "--gr") {
				Configurator::getInstance().gen_robust = true;
				cout << "Will also create robust model." << endl;
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

			if (input == "--logy") {
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

			if (input == "--min_log_range" || input == "--milog") {
				if (argc <= i) {
					std::cerr << "Missing argument for parameter min_log_range. Terminating." << std::endl;
					exit(-1);
				}
				Configurator::getInstance().min_log_range = atof(argv[i + 1]);
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
					std::cerr << "Missing argument for parameter min_pol_range. Terminating." << std::endl;
					exit(-1);
				}

				Configurator::getInstance().min_pol_range = atof(argv[i + 1]);
				i++;
			}

			if (input == "--min_log_range" || input == "--milog") {
				if (argc <= i) {
					std::cerr << "Missing argument for parameter min_log_range. Terminating." << std::endl;
					exit(-1);
				}
				cout << string(argv[i + 1]) << endl;
				Configurator::getInstance().min_log_range = atof(argv[i + 1]);
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


			// Factorial model
			if (input == "--fac_pol_min_coeff" || input == "--fpmic") {
				if (argc <= i) {
					std::cerr << "Missing argument for parameter fac_pol_min_coeff. Terminating." << std::endl;
					exit(-1);
				}
				Configurator::getInstance().min_fac_pol_range = atof(argv[i + 1]);
				i++;
			}

			if (input == "--fac_pol_max_coeff" || input == "--fpmac") {
				if (argc <= i) {
					std::cerr << "Missing argument for parameter fac_pol_max_coeff. Terminating." << std::endl;
					exit(-1);
				}
				Configurator::getInstance().max_fac_pol_range = atof(argv[i + 1]);
				i++;
			}

			if (input == "--fac_log_min_coeff" || input == "--flmic") {
				if (argc <= i) {
					std::cerr << "Missing argument for parameter fac_log_min_coeff. Terminating." << std::endl;
					exit(-1);
				}
				Configurator::getInstance().min_fac_log_range = atof(argv[i + 1]);
				i++;
			}

			if (input == "--fac_log_max_coeff" || input == "--flmac") {
				if (argc <= i) {
					std::cerr << "Missing argument for parameter fac_log_max_coeff. Terminating." << std::endl;
					exit(-1);
				}
				Configurator::getInstance().max_fac_log_range = atof(argv[i + 1]);
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

		if (Configurator::getInstance().texfile == "")
			Configurator::getInstance().texfile = "DefaultModel";

		if (Configurator::getInstance().outpath == "")
			Configurator::getInstance().outpath = "C:\\temp";

		return 0;
	}

	int SimAnMo::readInputMeasurementData(MeasurementDB*& inputDB) {
		std::string inputfile = Configurator::getInstance().inputfile;
		MeasurementDBReader dbreader = MeasurementDBReader();
		try {
			inputDB = dbreader.readInputFile(inputfile);
		}
		catch (...) {
			return 1;
		}
		return 0;
	}

//} // End Namespace
