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
	double max_pol_range;
	double max_log_range;
	std::string inputfile;
	std::string texfile;
	std::string outpath;

private:
	Configurator() {
		texfile = "";
		outpath = "";
		std_exp_range = 20;
		num_threads = 1;
		max_pol_range = 3;
		max_log_range = 2.5;
	}

	Configurator(const Configurator&);
	~Configurator() {}



};

#endif // !SOLUTIONCONFIGURATOR_H

