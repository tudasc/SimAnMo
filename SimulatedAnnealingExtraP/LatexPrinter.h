#ifndef LATEXPRINTER_H
#define LATEXPRINTER_H

#include "Solution.h"
#include "MeasurementDB.h"
#include <cstring>
#include <string.h>
#include <fstream>
#include "LinearSolution.h"

using namespace std;

template<class SolType>
struct CalcuationInfo{
	AbstractSolution* min_sol = NULL;
	double runtime = 0.0;
	double RSScost = 0.0;
	int thread_with_solution = -1;
	unsigned int iterations = 0;
	std::vector<SolType> sol_per_thread;
	MeasurementDB* datapoints;
	LinearSolution lin_sol;
	bool print_lin_sol = false;
	SolType ref_solution;
	bool print_ref_solution = false;
	bool print_measurepoints = false;
};

template<class SolType>
class LatexPrinter {
public: 
	LatexPrinter() {}
	void printSolution(std::string filename, AbstractSolution* sol, MeasurementDB* mdb, CalcuationInfo<SolType>& calcinf);
	void printColorDefinitions(int number, ofstream&  stream );
private:
	string colStr(int no);
	void printPrediction(ofstream&, CalcuationInfo<SolType>&, int stepsize);
	void printCostDevelopment(ofstream&, int);
	void printCostDetails(ofstream&, CalcuationInfo<SolType>&, int stepsize);
};

#endif