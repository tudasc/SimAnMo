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
	inline void noLogModel() { this->create_log_exp_model = false; }

	// Configuration variables
	// For Standard Solution
	double std_exp_range;
	int num_threads;

	double min_pol_range;
	double max_pol_range;

	double min_log_range;
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
		std_exp_range = 3.1;
		num_threads = 1;

		min_pol_range = 0;
		max_pol_range = 8;

		min_log_range = 0.0;
		max_log_range = 2.5;

		create_log_exp_model = true;
		print_confidence = false;
		confidence_interval = 0.0;
	}

	Configurator(const Configurator&);
	~Configurator() {}



};

#endif // !SOLUTIONCONFIGURATOR_H

