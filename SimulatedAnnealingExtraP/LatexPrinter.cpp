#include "LatexPrinter.h"
#include <iostream>
#include <fstream>
#include "QualityLogger.h"
#include <string> 
#include "ExtraPSolution.h"
#include "ExponentialSolution.h"
#include "ExponentialPolynomSolution.h"
#include "FactorialSolution.h"
#include "Configurator.h"
#include <cmath>

using namespace std;

template<class SolType>
void LatexPrinter<SolType>::printSolution(std::string filename, AbstractSolution* sol, MeasurementDB * mdb, CalcuationInfo<SolType>& calcinf)
{
	string rm_str = "";
	string slash_str = "";
#ifdef	__linux__
	rm_str = "rm";
	slash_str = "/";
	//cout << "Using Linux" << endl;
	//cout << Configurator::getInstance().outpath << endl;
#else
	//cout << "Using Windows" << endl;
	rm_str = "del";
	slash_str = "\\\\";
#endif

	calcinf.min_sol = sol;
	ofstream myfile;
	std::string filepath = Configurator::getInstance().outpath + slash_str + Configurator::getInstance().texfile + ".tex";
	myfile.open(filepath);

	myfile << "\\documentclass{article}" << endl
		<< "\\usepackage[letterpaper, margin=1.5cm]{geometry}" << endl
		<< "\\usepackage{xcolor, colortbl}" << endl
		<< "\\usepackage{pgfplots}" << endl
		<< "\\pgfplotsset{compat=newest}" << endl
		<< "\\usetikzlibrary{pgfplots.statistics}" << endl
		<< "\\usetikzlibrary{patterns}" << endl
		<< "\\usetikzlibrary{calc}" << endl
		<< "\\usepgfplotslibrary{fillbetween}" << endl;

	printColorDefinitions(Configurator::getInstance().num_threads, myfile);


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
		<< "\\item Cost ( " << sol->_cost_calc_type
		<< " ):" << calcinf.RSScost << endl
		<< "\\item Metrics of our model: " << sol->metricDetailsToStr();

	if (Configurator::getInstance().create_log_exp_model) {
		myfile << "\\item Metrics of lin-log model : " << calcinf.min_sol_log->metricDetailsToStr();
	}

	myfile << "\\end{itemize}" << endl
			<< "}" << endl;

	// Print Information on reference model if configured
	if (calcinf.print_ref_solution && 1 == 2) {
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
		<< "xlabel={rank of lattice $n$}," << endl
		<< "legend style={at={(0.03, 0.999)},anchor=north west, legend cell align=left, font=\\footnotesize, fill=none, draw=none}," << endl
		<< "]" << endl;

	// print the found model function
	myfile << "\\addplot [domain= " << mdb->getPairAt(0).first << ":" << mdb->getPairAt(mdb->get_size() - 1).first
		<< ", samples=110,unbounded coords=jump, draw=blue, very thick] " << endl
		<< sol->printModelFunctionLatex().c_str() << ";" << endl
		<< "\\addlegendentry{Our model $"
		<<	sol->printModelFunctionLatexShow()
		<< "$ };" << endl;

	// print the logarithmized model function if configured
	if (Configurator::getInstance().create_log_exp_model)
	{
		myfile << "\\addplot [draw=darkorange , domain= " << mdb->getPairAt(0).first << ":" << mdb->getPairAt(mdb->get_size() - 1).first
			<< ", samples=110,unbounded coords=jump, very thick] " << endl
			<< calcinf.min_sol_log->printModelFunctionLatex(0.0, true).c_str() << ";" << endl
			<< "\\addlegendentry{Lin-Log model $"
			<< calcinf.min_sol_log->printModelFunctionLatexShow(true)
			<< "$ };" << endl;
	}

	// print the reference model function if configured
	if (calcinf.print_ref_solution && 1==2) {
		myfile << "\\addplot [domain=1:" << mdb->getPairAt(mdb->get_size() - 1).first
			<< ", samples=110,unbounded coords=jump, draw=magenta, very thick] " << endl
			<< calcinf.ref_solution.printModelFunctionLatex().c_str() << ";" << endl
			<< "\\addlegendentry{Reference model};" << endl;
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

	if (Configurator::getInstance().create_log_exp_model)
	{
		printLogModel(myfile, sol, mdb, calcinf);
	}

// If configured, also show the Exponential Model
// Print powed logarithmic model
/*	if (Configurator::getInstance().create_log_exp_model) {
		myfile << "\\addplot [name path=func, domain=" << xstart << ":" << xdimension
			<< ", samples=110,unbounded coords=jump, draw=color6, very thick] " << endl
			<< cinfo.min_sol_log->printModelFunctionLatex(0.0, true).c_str() << ";" << endl
			<< "\\addlegendentry{Model from Log2};" << endl;
	}*/

/////////////////// Print the comparison of our and linear model
	int stepsize = QualityLogger::getInstance().get_size() / 200;
	printPrediction(myfile, calcinf, stepsize, sol);


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
		<< "ymode = log," << endl
		<< "ylabel={costs}," << endl
		<< "xlabel={step}," << endl
		<< "legend style={at={(0.99, 0.999)},anchor=north east, legend cell align=left, font=\\footnotesize, fill=none, draw=none}," << endl
		<< "]" << endl;

	if(Configurator::getInstance().print_costs)
		printCostDevelopment(myfile, stepsize);

	myfile << "\\end{axis}" << endl;
	myfile << "\\end{tikzpicture}" << endl
		<< "\\end{figure}" << endl;

//////////////////// Print more details
	if (Configurator::getInstance().print_costs)
		printCostDetails(myfile, calcinf, stepsize);

	myfile	<< "\\end{document}" << endl;


	myfile.close();

	// LaTeX Config: Adapt commands 
	std::string close_command = "cd " + Configurator::getInstance().path_pdf_xchange +
		" && \"PDFXCview.exe\" /close \"" + Configurator::getInstance().outpath + slash_str + Configurator::getInstance().texfile + ".pdf\"";	
#ifdef _WIN32
	std::cout << close_command << std::endl;
	system(close_command.c_str());
#endif

	std::string create_command = "cd " + Configurator::getInstance().outpath + " && pdflatex " + Configurator::getInstance().texfile + ".tex"; //   -interaction=batchmode
	cout << create_command << endl;
	system(create_command.c_str());
	std::string open_command = "cd " + Configurator::getInstance().path_pdf_xchange + 
		" && \"PDFXCview.exe\" /A \"page=1&zoom=125\" \"" + Configurator::getInstance().outpath + slash_str + Configurator::getInstance().texfile + ".pdf\"";
	
	std::string clean_command = "cd " + Configurator::getInstance().outpath + " && " + rm_str + " " + Configurator::getInstance().texfile + ".log " + Configurator::getInstance().texfile + ".aux";
	system(clean_command.c_str());	


#ifdef _WIN32
	std::cout << clean_command << std::endl;
	std::cout << open_command << std::endl;	
	if(Configurator::getInstance().open_latex_output)
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

	for (int i = 8; i < Configurator::getInstance().num_threads; i++) {
		colors.push_back("\\definecolor{color" + std::to_string(i+1) +
			"}{RGB}{" + std::to_string(rand() % 255) + ","
			+ std::to_string(rand() % 255) + ","
			+ std::to_string(rand() % 255) + "}");
	}

	for (size_t i = 0; i < colors.size(); i++) {
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
void LatexPrinter<SolType>::printLogModel(ofstream & myfile, AbstractSolution* sol, MeasurementDB * mdb, CalcuationInfo<SolType>& calcinf)
{
	int base = Configurator::getInstance().base_for_lin_log;
	myfile << endl << "\\newpage" << endl
		<< "\\section{Log-Function Evaluation}" << endl;

	myfile << "{\\large" << endl
		<< "\\begin{equation}" << endl
		<< "f_{log-model}(x) = "
		<< calcinf.min_sol_log->printModelFunctionLatexShow().c_str() << std::endl
		<< "\\end{equation}}" << endl;

	myfile << "\\begin{figure}[htb]" << endl
		<< "\\centering" << endl
		<< "\\setlength\\figureheight{10cm}" << endl
		<< "\\setlength\\figurewidth{10cm}" << endl
		<< "\\pgfplotsset{every tick label/.append style={font=\\small}}" << endl
		<< "\\begin{tikzpicture}" << endl
		<< "\\begin{axis}[" << endl
		<< "width=\\figurewidth," << endl
		<< "height=\\figureheight," << endl
		<< "scale only axis," << endl
		<< "ylabel={runtime in s}," << endl
		<< "xlabel={rank of lattice $n$}," << endl
		<< "legend style={at={(0.03, 0.999)},anchor=north west, legend cell align=left, font=\\footnotesize, fill=none, draw=none}," << endl
		<< "]" << endl;

	// print the found model function
	myfile << "\\addplot [domain= " << mdb->getPairAt(0).first << ":" << mdb->getPairAt(mdb->get_size() - 1).first
		<< ", samples=110,unbounded coords=jump, draw=darkorange, very thick] " << endl
		<< calcinf.min_sol_log->printModelFunctionLatex().c_str() << ";" << endl
		<< "\\addlegendentry{Lin-log model $"
		<< calcinf.min_sol_log->printModelFunctionLatexShow()
		<< "$ };" << endl;

	// Draw all logarithmized training data points
	myfile << "\\addplot[only marks, mark=diamond*,mark options={scale=1.5, fill=red},draw=red] coordinates {	" << endl;
	for (int i = 0; i < mdb->get_size(); i++) {
		pair<double, double> act_pair = mdb->getPairAt(i);
		myfile << "(" << act_pair.first << ", " << log(act_pair.second)/log(base) << ")" << endl;
	}
	myfile << "};" << endl
		<< "\\addlegendentry{ Training data };" << endl;

	myfile << "\\end{axis}" << endl
		<< "\\end{tikzpicture}" << endl
		<< "\\end{figure}" << endl;


}


template<class SolType>
void LatexPrinter<SolType>::printPrediction(ofstream & myfile, CalcuationInfo<SolType> & cinfo, int stepsize, AbstractSolution* sol)
{
	int no_points = cinfo.datapoints->get_size();
	AbstractSolution& minsol = cinfo.sol_per_thread[cinfo.thread_with_solution];
	double cv = Configurator::getInstance().confidence_interval;

	int xdimension = (int)((cinfo.datapoints->getPairAt(no_points - 1).first) * 1.00);
	int xstart = (int)(cinfo.datapoints->getPairAt(0).first);

	if (cinfo.datapoints->get_measures_size() > 0)
	{
		std::pair<double, double> act_pair = cinfo.datapoints->getMeasurePairAt(0);
		xstart = std::max<int>((int)act_pair.first-6, 0);
	}


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
				cinfo.datapoints->getMeasurePairAt(cinfo.datapoints->get_measures_size() - 1).first);
		}
	}

	double funcVal = sol->evaluateModelFunctionAt(xdimension);
	double funcValMinus = sol->evaluateModelFunctionAt(xdimension, -cv);
	double maxaddedpoint = -1;
	if (cinfo.datapoints->get_measures_size() > 0) {
		maxaddedpoint = cinfo.datapoints->getMeasurePairAt(cinfo.datapoints->get_measures_size() - 1).second;
	}
	funcVal = std::max<double>(funcVal, maxaddedpoint);
	
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
		<< "scale only axis," << endl;
	if (Configurator::getInstance().ymode_log == true) {
		myfile << "ymode=log," << endl;
	}		
	myfile << "xlabel={rank of lattice $n$}," << endl
		<< "ylabel={runtime in s}," << endl
		<< "legend style={at={(0.03, 0.999)},anchor=north west, legend cell align=left, font=\\footnotesize, fill=none, draw=none}," << endl
		<< "]" << endl;

	// print the found linear function
	/*myfile << "\\addplot [domain=1:" << xdimension
		<< ", samples=110,unbounded coords=jump, draw=green, very thick] " << endl;
	myfile << "(\\x, {2 ^ ( " << cinfo.lin_sol.getAt(0)
		<< " + " << cinfo.lin_sol.getAt(1) << " * \\x )) });" << endl
		<< "\\addlegendentry{Linear Model};" << endl;*/

	// print the found model function
	myfile << "\\addplot [name path=func, domain=" << xstart << ":" << xdimension
		<< ", samples=110,unbounded coords=jump, draw=blue, very thick] " << endl
		<< minsol.printModelFunctionLatex().c_str() << ";" << endl
		<< "\\addlegendentry{Our model $"
		<< minsol.printModelFunctionLatexShow()
		<< "$ };" << endl;

	if (Configurator::getInstance().print_confidence) {		
		myfile << "\\addplot [name path=funcplus, domain=" << xstart << ":" << xdimension
			<< ", samples=110,unbounded coords=jump, draw=color5, dashed, ] " << endl
			<< minsol.printModelFunctionLatex(cv).c_str() << ";" << endl
			<< "\\addlegendentry{Model + " << cv << "};" << endl;

		myfile << "\\addplot [name path=funcminus, domain=" << xstart << ":" << xdimension
			<< ", samples=110,unbounded coords=jump, draw=red, dashed, ] " << endl
			<< minsol.printModelFunctionLatex(-cv).c_str() << ";" << endl
			<< "\\addlegendentry{Model - " << cv << "};" << endl;
	}

	// print the reference solution if configured
	if(cinfo.print_ref_solution && 1==2) {
		myfile << "\\addplot [domain=" << xstart << ":" << xdimension
			<< ", samples=110, unbounded coords=jump, draw=magenta, very thick] " << endl
			<< cinfo.ref_solution.printModelFunctionLatex().c_str() << ";" << endl
			<< "\\addlegendentry{Reference model};" << endl;
	}

	// Print powed logarithmic model if configured
	if (Configurator::getInstance().create_log_exp_model) {
		myfile << "\\addplot [name path=func, domain=" << xstart << ":" << xdimension
			<< ", samples=110,unbounded coords=jump, draw=darkorange, very thick] " << endl
			<< cinfo.min_sol_log->printModelFunctionLatex(0.0, true).c_str() << ";" << endl
			<< "\\addlegendentry{Lin-log model $"
			<< cinfo.min_sol_log->printModelFunctionLatexShow(true)
			<< "$ };" << endl;
	}

	// print additional measurement points if configured
	if (cinfo.print_measurepoints) {
		myfile << "\\addplot[only marks, mark = *, mark options = {scale=1.5, fill=color5}, draw = color5] "
			<< "coordinates { " << std::endl;
		for (int i = 0; i < cinfo.datapoints->get_measures_size(); i++) {
			std::pair<double, double> act_pair = cinfo.datapoints->getMeasurePairAt(i);
			myfile << "( " << act_pair.first << " , " << act_pair.second << " )" << std::endl;
		}
		myfile << "};" << std::endl
			<< "\\addlegendentry{ Measurements };" << endl;

		/*for (int i = 0; i < cinfo.datapoints->get_measures_size(); i++) {
			std::pair<double, double> act_pair = cinfo.datapoints->getMeasurePairAt(i);
			double start_y = act_pair.second - (funcVal*0.43);
			double end_y = act_pair.second + (funcVal*0.43);
			myfile << "\\addplot +[mark=none, thick, black, dashed] coordinates {(" << act_pair.first 
				<< ", " << start_y << ") ( " << act_pair.first << ", " << end_y << ")};" << endl;
		}*/
		if (Configurator::getInstance().print_confidence) {
			myfile << "\\draw[help lines, name path = clippath, white, very thin]"
				<< "(" << xstart << ", " << funcValMinus << ") -- " << "(" << xdimension << " , " << funcValMinus << ");" << endl;

			myfile << "\\addplot fill between[" << endl
				<< "of = funcplus and funcminus," << endl
				<< "every even segment/.style = { gray,opacity = .2 }," << endl
				<< "soft clip = { clippath }," << endl
				<< "];" << endl;

			myfile << "\\addplot fill between[" << endl
				<< "of = funcplus and funcminus," << endl
				<< "every even segment/.style = { gray,opacity = .2 }," << endl
				<< "soft clip = { domain = " << xstart << ":" << xdimension << "}," << endl
				<< "]; " << endl;
		}
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
			<< " is " << ourdef << "\\%. (" << oury << " vs. " << measurey << ")" << std::endl;

		
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
		for (int j = 1; j < QualityLogger::getInstance().get_size(tid); j += stepsize) {
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

	for (int tid = 0; tid < (int)calcinf.sol_per_thread.size(); tid++) {
		const SolType & sol = calcinf.sol_per_thread[tid];
		myfile << "\\subsection{Thread " << tid << "}" << endl
			<< "\\begin{itemize}" << endl
			<< "\\item Cost ( " << sol._cost_calc_type << " ):" << calcinf.sol_per_thread[tid].get_costs() << endl
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
			<< "ymode = log," << endl
			<< "ymax=" << QualityLogger::getInstance().get_max_cost_to_print(stepsize, tid) * 1.1 << "," << endl
			<< "ylabel={costs}," << endl
			<< "xlabel={step}," << endl
			<< "legend style={at={(0.99, 0.999)},anchor=north east, legend cell align=left, font=\\footnotesize, fill=none, draw=none}," << endl
			<< "]" << endl;

			myfile << "\\addplot[thick] coordinates {	" << endl;

			// Run over all Log-Entries and print them
			for (int j = 1; j < QualityLogger::getInstance().get_size(tid); j += stepsize) {
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
template class LatexPrinter<ExponentialSolution>;
template class LatexPrinter<ExponentialPolynomSolution>;
template class LatexPrinter<FactorialSolution>;
