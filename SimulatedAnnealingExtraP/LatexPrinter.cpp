#include "LatexPrinter.h"
#include <iostream>
#include <fstream>
#include "QualityLogger.h"
#include <string> 
#include "ExtraPSolution.h"
#include "Configurator.h"

using namespace std;

template<class SolType>
void LatexPrinter<SolType>::printSolution(std::string filename, AbstractSolution* sol, MeasurementDB * mdb, CalcuationInfo<SolType>& calcinf)
{
	calcinf.min_sol = sol;
	ofstream myfile;
	std::string filepath = Configurator::getInstance().outpath + "\\\\" + Configurator::getInstance().texfile + ".tex";
	myfile.open(filepath);

	myfile << "\\documentclass{article}" << endl
		<< "\\usepackage[letterpaper, margin=1.5cm]{geometry}" << endl
		<< "\\usepackage{xcolor, colortbl}" << endl
		<< "\\usepackage{pgfplots}" << endl
		<< "\\pgfplotsset{compat=newest}" << endl
		<< "\\usetikzlibrary{pgfplots.statistics}" << endl
		<< "\\usetikzlibrary{patterns}" << endl
		<< "\\usetikzlibrary{calc}" << endl;

	printColorDefinitions(8, myfile);


	myfile << "\\begin{document}" << endl
		<< "\\section{Evaluation}" << endl;

	// Print Information on found model
	myfile
		<< "{\\large" << endl
		<< "\\begin{equation}" << endl
		<< "f_{model}(x) = "
		<< sol->printModelFunctionLatexShow().c_str() << std::endl
		<< "\\end{equation}"
		<< "\\begin{itemize}" << endl
		<< "\\item Thread " << calcinf.thread_with_solution << " found solution in "
		<< calcinf.iterations << " steps requiring "
		<< calcinf.runtime << " seconds." << endl
		<< "\\item Cost (RSS):" << calcinf.RSScost << endl
		<< "\\end{itemize}" << endl;
	myfile << "}" << endl;

	// Print Information on reference model if configured
	if (calcinf.print_ref_solution && 1==2) {
		myfile
		<< "\\begin{equation}" << endl
			<< "f_{model}(x) = "
			<< calcinf.ref_solution.printModelFunctionLatexShow().c_str() << std::endl
			<< "\\end{equation}"
			<< "\\begin{itemize}" << endl
			<< "\\item Cost (RSS):" << std::to_string(calcinf.ref_solution.get_costs()) << endl
			<< "\\end{itemize}" << endl;
	}

	/*myfile << "$f_{linear}(x)" << calcinf.ref_solution.getAt(0)
		<< " + " << calcinf.ref_solution.getAt(1) << " * 2 ^ {"
		<< calcinf.ref_solution.getAt(2) << "* x ^ {" << calcinf.ref_solution.getAt(3) << "}}$" << endl;*/

	/*if (calcinf.print_lin_sol) {
	myfile 	<< "\\begin{equation}" << endl
		<< "f_{linear}(x) = 2 ^ {" << calcinf.lin_sol.getAt(0)
		<< " + " << calcinf.lin_sol.getAt(1) << " * x}\\\\" << endl
		<< "\\end{equation}";
	}*/		
	myfile << "\\begin{figure}[htb]" << endl
		<< "\\centering" << endl
		<< "\\newlength\\figureheight" << endl
		<< "\\newlength\\figurewidth" << endl
		<< "\\setlength\\figureheight{10cm}" << endl
		<< "\\setlength\\figurewidth{10cm}" << endl
		<< "\\pgfplotsset{every tick label/.append style={font=\\small}}" << endl
		<< "\\begin{tikzpicture}" << endl
		<< "\\begin{axis}[" << endl
		<< "width=\\figurewidth," << endl
		<< "height=\\figureheight," << endl
		<< "scale only axis," << endl
		//<< "ymax=" << mdb->getPairAt(mdb->get_size() - 1).second << "," << endl
		<< "ylabel={runtime in s}," << endl
		<< "xlabel={rank of lattice d}," << endl
		<< "legend style={at={(0.20, 0.999)},anchor=north west, legend cell align=left}," << endl
		<< "]" << endl;

	// print the found model function
	myfile  << "\\addplot [domain= " << mdb->getPairAt(0).first << ":" <<  mdb->getPairAt(mdb->get_size()-1).first
		<< ", samples=110,unbounded coords=jump, draw=blue, very thick] " << endl
		<< sol->printModelFunctionLatex().c_str() << ";" << endl
		<< "\\addlegendentry{Our Model};" << endl;

	// print the reference model function if configured
	if (calcinf.print_ref_solution && 1==2) {
		myfile << "\\addplot [domain=1:" << mdb->getPairAt(mdb->get_size() - 1).first
			<< ", samples=110,unbounded coords=jump, draw=magenta, very thick] " << endl
			<< calcinf.ref_solution.printModelFunctionLatex().c_str() << ";" << endl
			<< "\\addlegendentry{Reference Model};" << endl;
	}


	// print the found linear function
	/*myfile << "\\addplot [domain=1:" << mdb->getPairAt(mdb->get_size() - 1).first
		<< ", samples=110,unbounded coords=jump, draw=green, very thick] " << endl;
	myfile << "(\\x, {2 ^ ( " << calcinf.lin_sol.getAt(0)
		<< " + " << calcinf.lin_sol.getAt(1) << " * \\x ) });" << endl
		<< "\\addlegendentry{Linear Model};" << endl;*/

	myfile	<< "\\addplot[only marks, mark=diamond*,mark options={scale=1.5, fill=red},draw=red] coordinates {	" << endl;
		// Run over all Measurements and print them
	for (int i = 0; i < mdb->get_size(); i++) {
		pair<double, double> act_pair = mdb->getPairAt(i);
		myfile << "(" << act_pair.first << ", " << act_pair.second << ")" << endl;
	}
	myfile << "};" << endl
		<< "\\addlegendentry{ Training data };" << endl;

	myfile << "\\end{axis}" << endl
		<< "\\end{tikzpicture}" << endl
		<< "\\end{figure}" << endl;

/////////////////// Print the comparison of our and linear model
	int stepsize = QualityLogger::getInstance().get_size() / 200;
	printPrediction(myfile, calcinf, stepsize);


/////////////////// Draw the quality development


	if (stepsize == 0)
		stepsize = 1;

	myfile << "\\newpage" << endl
		<< "\\section{Development of Costs}" << endl
		<< "\\begin{figure}[htb]" << endl
		<< "\\centering" << endl
		<< "\\setlength\\figureheight{10cm}" << endl
		<< "\\setlength\\figurewidth{10cm}" << endl
		<< "\\begin{tikzpicture}" << endl
		<< "\\begin{axis}[" << endl
		<< "width=\\figurewidth," << endl
		<< "height=\\figureheight," << endl
		<< "scale only axis," << endl
		//<< "ymax=" << QualityLogger::getInstance().get_max_cost() * 1.2 << "," << endl
		//<< "ymax=" << QualityLogger::getInstance().get_max_cost_to_print(stepsize) * 1.1 << "," << endl
		<< "ylabel={costs}," << endl
		<< "xlabel={step}," << endl
		<< "legend style={at={(0.99, 0.999)},anchor=north east, legend cell align=left}," << endl
		<< "]" << endl;

	printCostDevelopment(myfile, stepsize);

	myfile << "\\end{axis}" << endl;
	myfile << "\\end{tikzpicture}" << endl
		<< "\\end{figure}" << endl;

//////////////////// Print more details
	printCostDetails(myfile, calcinf, stepsize);

	myfile	<< "\\end{document}" << endl;


	myfile.close();

	string rm_str;
#ifdef	__linux__
	rm_str = "rm";
#else
	rm_str = "del";
#endif

	// LaTeX Config: Adapt commands 
	std::string close_command = "cd C:\\Program Files\\Tracker Software\\PDF Viewer && \"PDFXCview.exe\" /close \"" + Configurator::getInstance().outpath + "\\\\" + Configurator::getInstance().texfile + ".pdf\"";
	std::cout << close_command << std::endl;
#ifdef _WIN32
	system(close_command.c_str());
#endif

	std::string create_command = "cd " + Configurator::getInstance().outpath + " && pdflatex " + Configurator::getInstance().texfile + ".tex"; //   -interaction=batchmode
	system(create_command.c_str());
	std::string open_command = "cd C:\\Program Files\\Tracker Software\\PDF Viewer && \"PDFXCview.exe\" /A \"page=1&zoom=125\" \"" + Configurator::getInstance().outpath + "\\\\"+ Configurator::getInstance().texfile + ".pdf\"";
	
	std::string clean_command = "cd " + Configurator::getInstance().outpath + " && " + rm_str + " " + Configurator::getInstance().texfile + ".log " + Configurator::getInstance().texfile + ".aux";
	system(clean_command.c_str());	
	std::cout << create_command << std::endl;
	std::cout << open_command << std::endl;	

#ifdef _WIN32
	system(open_command.c_str());
#endif
	system("exit");

	exit( 0 );
}

template<class SolType>
void LatexPrinter<SolType>::printColorDefinitions(int number, ofstream & stream)
{
	vector<string> colors;
	colors.push_back("\\definecolor{color1}{RGB}{255,0,0}");
	colors.push_back("\\definecolor{color2}{RGB}{0,255,0}");
	colors.push_back("\\definecolor{color3}{RGB}{0,0,255}");
	colors.push_back("\\definecolor{color4}{RGB}{255,255,0}");
	colors.push_back("\\definecolor{color5}{RGB}{255,0,255}");
	colors.push_back("\\definecolor{color6}{RGB}{0,255,255}");
	colors.push_back("\\definecolor{color7}{RGB}{90,90,255}");
	colors.push_back("\\definecolor{color8}{RGB}{90,255,90}");

	for (int i = 0; i < number; i++) {
		stream << colors[i].c_str() << std::endl;
	}	

	stream << "\\definecolor{darkorange}{RGB}{255,140,0}" << std::endl;
}

template<class SolType>
string LatexPrinter<SolType>::colStr(int no)
{
	string str = "color" + std::to_string(no);
	return str;
}

template<class SolType>
void LatexPrinter<SolType>::printPrediction(ofstream & myfile, CalcuationInfo<SolType> & cinfo, int stepsize)
{
	int no_points = cinfo.datapoints->get_size();
	AbstractSolution& minsol = cinfo.sol_per_thread[cinfo.thread_with_solution];

	int xdimension = (int)((cinfo.datapoints->getPairAt(no_points - 1).first) * 1.03);

	double ourfuncmax = minsol.evaluateModelFunctionAt(xdimension);
	double linearfuncmax = cinfo.lin_sol.evaluateModelFunctionAt(xdimension);

	double ydimension = ourfuncmax;

	if (cinfo.print_lin_sol) {
		ydimension = std::max<double>
			(ourfuncmax, linearfuncmax);
	}

	if (cinfo.print_measurepoints) {
		if (cinfo.datapoints->get_measures_size() > 0) {
			ydimension = (int)std::max<double>(ydimension,
				cinfo.datapoints->getMeasurePairAt(cinfo.datapoints->get_measures_size() - 1).second);
			xdimension = (int)std::max<double>(xdimension,
				cinfo.datapoints->getMeasurePairAt(cinfo.datapoints->get_measures_size() - 1).first) + 1;
		}
	}

	myfile << "\\newpage" << endl;
	myfile << "\\section{Prediction of Runtime}" << endl;

	myfile << "\\begin{figure}[htb]" << endl
		<< "\\centering" << endl
		<< "\\setlength\\figureheight{10cm}" << endl
		<< "\\setlength\\figurewidth{10cm}" << endl
		<< "\\begin{tikzpicture}" << endl
		<< "\\begin{axis}[" << endl
		<< "width=\\figurewidth," << endl
		<< "height=\\figureheight," << endl
		<< "scale only axis," << endl
		//<< "ymax=" << ydimension * 1.05 << ","
		<< "ylabel={costs}," << endl
		<< "xlabel={step}," << endl
		<< "legend style={at={(0.05, 0.95)},anchor=north west, legend cell align=left}," << endl
		<< "]" << endl;

	// print the found linear function
	/*myfile << "\\addplot [domain=1:" << xdimension
		<< ", samples=110,unbounded coords=jump, draw=green, very thick] " << endl;
	myfile << "(\\x, {2 ^ ( " << cinfo.lin_sol.getAt(0)
		<< " + " << cinfo.lin_sol.getAt(1) << " * \\x )) });" << endl
		<< "\\addlegendentry{Linear Model};" << endl;*/

	// print the found model function
	myfile << "\\addplot [domain=1:" << xdimension
		<< ", samples=110,unbounded coords=jump, draw=blue, very thick] " << endl
		<< minsol.printModelFunctionLatex().c_str() << ";" << endl
		<< "\\addlegendentry{Our Model};" << endl;

	// print the reference solution if configured
	if(cinfo.print_ref_solution && 1==2) {
		myfile << "\\addplot [domain=1:" << xdimension
			<< ", samples=110, unbounded coords=jump, draw=magenta, very thick] " << endl
			<< cinfo.ref_solution.printModelFunctionLatex().c_str() << ";" << endl
			<< "\\addlegendentry{Reference Model};" << endl;
	}

	// print additional measurement points if configured
	if (cinfo.print_measurepoints) {
		myfile << "\\addplot[only marks, mark = *, mark options = {scale=1.5, fill=darkorange}, draw = darkorange] "
			<< "coordinates { " << std::endl;
		for (int i = 0; i < cinfo.datapoints->get_measures_size(); i++) {
			std::pair<double, double> act_pair = cinfo.datapoints->getMeasurePairAt(i);
			myfile << "( " << act_pair.first << " , " << act_pair.second << " )" << std::endl;
		}
		myfile << "};" << std::endl
			<< "\\addlegendentry{ Measurements };";
	}

	myfile << "\\end{axis}" << endl
		<< "\\end{tikzpicture}" << endl
		<< "\\end{figure}" << endl;

	if (cinfo.print_ref_solution) {
		double measurex = (double)xdimension;
		double oury = cinfo.min_sol->evaluateModelFunctionAt(measurex);
		double refy = cinfo.ref_solution.evaluateModelFunctionAt(measurex);
		double ourrefdef = (abs(oury - refy) / (refy)) * 100.0;

		myfile << "{\\large"
			<< "\\begin{itemize}" << std::endl
			<< "\\item " << "Deviation of our model and reference at x=" << measurex
			<< " is " << ourrefdef << "\\%." << std::endl;
		myfile << "\\end{itemize}}" << std::endl;
	}

	// Print further information
	if (cinfo.datapoints->get_measures_size() > 0) {
		double measurex = cinfo.datapoints->getMeasurePairAt(cinfo.datapoints->get_measures_size() - 1).first;
		double measurey = cinfo.datapoints->getMeasurePairAt(cinfo.datapoints->get_measures_size() - 1).second;
		double oury = cinfo.min_sol->evaluateModelFunctionAt(measurex);
		double refy = cinfo.ref_solution.evaluateModelFunctionAt(measurex);
		double ourdef = (abs(oury - measurey) / (measurey)) * 100.0;

		myfile << "{\\large"
			<< "\\begin{itemize}" << std::endl
			<< "\\item " << "Deviation of our model and measures at x=" << measurex
			<< " is " << ourdef << "\\%." << std::endl;

		
		if (cinfo.print_ref_solution) {			
			double refdef = (abs(refy - measurey) / (measurey)) * 100.0;
			myfile << "\\item " << "Deviation of reference model and measures at x=" << measurex
				<< " is " << refdef << "\\%." << std::endl;
		}

		myfile << "\\end{itemize}}" << std::endl;
	}
}

template<class SolType>
void LatexPrinter<SolType>::printCostDevelopment(ofstream & myfile, int stepsize)
{
	for (int tid = 0; tid < QualityLogger::getInstance().get_no_of_vectors(); tid++) {
		myfile << "\\addplot[only marks, mark=diamond,mark options={}," << endl
			<< "draw=" << colStr(tid + 1).c_str() << ",fill=" << colStr(tid + 1).c_str() << "] coordinates {	" << endl;

		// Run over all Log-Entries and print them
		for (int j = 0; j < QualityLogger::getInstance().get_size(tid); j += stepsize) {
			pair<unsigned int, double> act_pair = QualityLogger::getInstance().get_entry_at(j, tid);
			myfile << "(" << act_pair.first << ", " << act_pair.second << ")" << endl;
		}

		string tmpstr = "Cost Thread " + std::to_string(tid);

		myfile << "};" << endl
			<< "\\addlegendentry{ " << tmpstr << "};" << endl;
	}		
}

template<class SolType>
void LatexPrinter<SolType>::printCostDetails(ofstream & myfile, CalcuationInfo<SolType>& calcinf, int stepsize)
{
	myfile << "\\newpage" << endl;
	myfile << "\\section{Details of Costs}" << endl;

	for (int tid = 0; tid < calcinf.sol_per_thread.size(); tid++) {
		const SolType & sol = calcinf.sol_per_thread[tid];
		myfile << "\\subsection{Thread " << tid << "}" << endl
			<< "\\begin{itemize}" << endl
			<< "\\item Cost (RSS):" << calcinf.sol_per_thread[tid].get_costs() << endl
			<< "\\item Thread " << tid << " found solution:"
			<< "{\\large" << endl;

		myfile << "\\begin{equation}" << endl
			<< sol.printModelFunctionLatexShow().c_str() << std::endl
			<< "\\end{equation}" << endl;

		/*myfile << "$"
			<< sol.printModelFunctionLatexShow().c_str() << "$" << std::endl;*/
		/*<< sol.getAt(0)
			<< " + " << sol.getAt(1) << " * 2 ^ {"
			<< sol.getAt(2) << "* x ^ {" << sol.getAt(3) << "}}$" << endl;*/

		myfile << "}" << endl // End large
			<< "\\end{itemize}" << endl;

		// Now print the diagram for the costs
		myfile << "\\begin{figure}[htb]" << endl
			<< "\\centering" << endl
			<< "\\setlength\\figureheight{8cm}" << endl
			<< "\\setlength\\figurewidth{8cm}" << endl
			<< "\\begin{tikzpicture}" << endl
			<< "\\begin{axis}[" << endl
			<< "width=\\figurewidth," << endl
			<< "height=\\figureheight," << endl
			<< "scale only axis," << endl			
			<< "ymax=" << QualityLogger::getInstance().get_max_cost_to_print(stepsize, tid) * 1.1 << "," << endl
			<< "ylabel={costs}," << endl
			<< "xlabel={step}," << endl
			<< "legend style={at={(0.99, 0.999)},anchor=north east, legend cell align=left}," << endl
			<< "]" << endl;

			myfile << "\\addplot[thick] coordinates {	" << endl;

			// Run over all Log-Entries and print them
			for (int j = 0; j < QualityLogger::getInstance().get_size(tid); j += stepsize) {
				pair<unsigned int, double> act_pair = QualityLogger::getInstance().get_entry_at(j, tid);
				myfile << "(" << act_pair.first << ", " << act_pair.second << ")" << endl;
			}

			string tmpstr = "Cost Thread " + std::to_string(tid);

			myfile << "};" << endl
				<< "\\addlegendentry{ " << tmpstr << "};" << endl;

		myfile << "\\end{axis}" << endl;
		myfile << "\\end{tikzpicture}" << endl
			<< "\\end{figure}" << endl;
		myfile << "\\newpage" << endl;
	}
}

template class LatexPrinter<Solution>;
template class LatexPrinter<ExtraPSolution>;
template class LatexPrinter<LinearSolution>;