#ifndef SOLUTIONCONFIGURATOR_H
#define SOLUTIONCONFIGURATOR_H
#include <cstring>

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
	int std_exp_range;
	int num_threads;

	double min_pol_range;
	double max_pol_range;

	double max_log_range;
	std::string inputfile;
	std::string texfile;
	std::string outpath;

	bool create_log_exp_model;

	// Print Configuration
	bool print_confidence;
	double confidence_interval;

private:
	Configurator() {
		texfile = "";
		outpath = "";
		std_exp_range = 20;
		num_threads = 1;

		min_pol_range = 1;
		max_pol_range = 16;

		max_log_range = 2.5;

		create_log_exp_model = false;
		print_confidence = false;
		confidence_interval = 0.0;
	}

	Configurator(const Configurator&);
	~Configurator() {}



};

#endif // !SOLUTIONCONFIGURATOR_H

