#ifndef SOLUTIONCONFIGURATOR_H
#define SOLUTIONCONFIGURATOR_H
#include <cstring>
#include <string>

//using namespace std;

class Configurator {
public:
	static Configurator& getInstance()
	{
		static Configurator instance;
		return instance;
	}

	static Configurator* _instance;

	// Configuration variables
	// For Standard Solution
	double std_exp_range;
	int num_threads;

	double min_pol_range;
	double max_pol_range;

	// exp-pol models
	double min_exp_coeff_range;
	double max_exp_coeff_range;

	double min_exp_exp_range;
	double max_exp_exp_range;

	double min_log_range;
	double max_log_range;
	std::string inputfile;
	std::string texfile;
	std::string outpath;

	int no_of_trials;

	bool create_log_exp_model;

	// Print Configuration
	bool do_latex_output;
	std::string path_pdf_xchange;
	bool open_latex_output;

	bool print_confidence;
	double confidence_interval;

	int base_for_lin_log;



private:
	Configurator() {
		texfile = "";
		outpath = "";
		std_exp_range = 1.8;
		num_threads = 1;

		min_pol_range = 0;
		max_pol_range = 8;

		min_log_range = 0.0;
		max_log_range = 2.5;

		// exp-pol models
		min_exp_coeff_range = 0.1;
		max_exp_coeff_range = 1.99;

		min_exp_exp_range = 0.1;
		max_exp_exp_range = 1.5;

		no_of_trials = 1;

		// Print Configuration
		do_latex_output = false;
		open_latex_output = false;
		path_pdf_xchange = "C:\\Program Files\\Tracker Software\\PDF Viewer";

		print_confidence = false;
		confidence_interval = 0.0;

		create_log_exp_model = false;
		base_for_lin_log=2;
	}

	Configurator(const Configurator&);
	~Configurator() {}



};

#endif // !SOLUTIONCONFIGURATOR_H

