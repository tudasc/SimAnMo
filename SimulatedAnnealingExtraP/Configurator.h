#ifndef SOLUTIONCONFIGURATOR_H
#define SOLUTIONCONFIGURATOR_H
#include <cstring>
#include <string>

//using namespace std;

#define TYPE_EIGENPARAMETER 0
#define TYPE_ALGLIBPARAMETER 1

#define TYPE_COST_RSS 0
#define TYPE_COST_RARSD 1

class Configurator {
public:
	static Configurator& getInstance()
	{
		static Configurator instance;
		return instance;
	}

	~Configurator() {}

	//static Configurator& _instance;

	/*unsigned long long getID() {
		return  this->glob_id++;
	}*/

	std::string solution_type;
	std::string costcalc_type;

	// Configuration variables
	// For Standard Solution
	double std_exp_range;
	int num_threads;

	// Configuration for annealing itself
	int ann_steps;
	double ann_cooling_rate;
	double ann_target_temp;
	int ann_steps_wo_mod;
	int ann_steps_backtrack;

	double min_pol_range;
	double max_pol_range;

	// exp models
	double max_exp_range;

	// exp-pol models
	double min_exp_coeff_range;
	double max_exp_coeff_range;

	double min_exp_exp_range;
	double max_exp_exp_range;

	// factorial models
	double min_fac_coeff_range;
	double max_fac_coeff_range;

	double min_log_range;
	double max_log_range;
	std::string inputfile;
	std::string texfile;
	std::string outpath;

	int no_of_trials;

	bool create_log_exp_model;
	bool gen_robust;

	// Print Configuration
	bool do_latex_output;
	std::string path_pdf_xchange;
	bool open_latex_output;
	bool ymode_log;

	bool print_confidence;
	double confidence_interval;

	int base_for_lin_log;

	bool print_costs;

	unsigned long long glob_id;

	double max_cost; // A cap for the costs within the Solultion initializer
	int param_est_typ;

private:
	Configurator() {
		solution_type = "extrapsolution";
		costcalc_type = "nnrrsscostcalculator";

		// Configuration for annealing itself
		ann_steps = 35;
		ann_cooling_rate = 0.998;
		ann_target_temp = 1e-14;
		ann_steps_wo_mod = 200000;
		ann_steps_backtrack = 200000;

		texfile = "";
		outpath = "";

		//std_exp_range = 1.8;
		
		num_threads = 1;

		// pol-log models
		min_pol_range = -0.5;
		max_pol_range = 5.00;

		min_log_range = -0.50;
		max_log_range = 1.50;

		// exp models
		max_exp_range = 3.0;

		// exp-pol models
		min_exp_coeff_range = 1e-3;
		max_exp_coeff_range = 2.5;

		min_exp_exp_range = 0.1;
		max_exp_exp_range = 2.5;

		// factorial models
		min_fac_coeff_range = 1e-3;
		max_fac_coeff_range = 1e-1;

		no_of_trials = 1;

		// Print Configuration
		do_latex_output = false;
		open_latex_output = false;
		path_pdf_xchange = "C:\\Program Files\\Tracker Software\\PDF Viewer";

		ymode_log = false;
		gen_robust = false;

		print_confidence = false;
		confidence_interval = 0.0;

		create_log_exp_model = false;
		base_for_lin_log=2;
		print_costs = false;


		glob_id = 1;

		max_cost = std::numeric_limits<double>::max();// 10e70;// std::numeric_limits<double>::max() * 10e-300;
		param_est_typ = TYPE_EIGENPARAMETER;
	}

	Configurator(const Configurator&);
	Configurator& operator = (const Configurator&) {
		return *this;
	};



};

#endif // !SOLUTIONCONFIGURATOR_H

