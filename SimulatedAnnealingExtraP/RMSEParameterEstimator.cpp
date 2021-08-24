#include "RMSEParameterEstimator.h"
#include "EigenParameterEstimator.h"
#include <iostream>

RMSEParameterEstimator::RMSEParameterEstimator()
{
	_mdb = nullptr;
}

RMSEParameterEstimator::RMSEParameterEstimator(MeasurementDB* mdb) : _mdb(mdb)
{

}

RMSEParameterEstimator::~RMSEParameterEstimator()
{
	_mdb = nullptr;
}

int RMSEParameterEstimator::estimateParameters(AbstractSolution* sol, double newrelerr) {
	int m = _mdb->get_size();
	int n = 2; // Nummber of variables

	RMSEParameterEstimatorAlglibParameter params;

	params.x_n.resize(m);
	params.y_n.resize(m);
	params.mdb = _mdb;
	params.sol = sol;

	// Fill the constant values in the vectors for 
	for (int i = 0; i < m; i++) {
		params.x_n[i] = _mdb->getPairAt(i).first;
		params.y_n[i] = _mdb->getPairAt(i).second;
	}

	EigenParameterEstimator eigpar = EigenParameterEstimator(_mdb);
	int ret = eigpar.estimateParameters(sol);
	if (ret > 0 && ret < 100) {
		cout << "General error " << endl;
		return ERR_GENERAL_INVALID;
	}

	double par1 = sol->getAt(0);
	double par2 = sol->getAt(1);

	real_1d_array v_init; v_init.setlength(2);
	v_init[0] = par1; v_init[1] = par2;
	real_1d_array s = "[1 , 1]"; // scale does not make difference
	real_2d_array c = "[[1, 0, 10e-12], [0, 1 , 10e-18]]";
	integer_1d_array ct = "[1,1]"; // >, >

	minbleicstate state;
	double epsg = 0;
	double epsf = 0;
	double epsx = 0.000000001;
	ae_int_t maxits = 0;

	minbleiccreate(v_init, state);
	minbleicsetlc(state, c, ct, 0);
	minbleicsetscale(state, s);
	minbleicsetcond(state, epsg, epsf, epsx, maxits);

	minbleicoptguardsmoothness(state);
	minbleicoptguardgradient(state, 0.00001);
	minbleicreport rep;

	try
	{
		alglib::minbleicoptimize(state, RMSEParameterEstimatorFunctionGradEval, nullptr, &params);
	}

	catch (...) {
		return ERR_GENERAL_INVALID;
	}

	minbleicresults(state, v_init, rep);
	//printf("%d\n", int(rep.terminationtype)); // EXPECTED: 4
	//printf("%s\n", v_init.tostring(2).c_str()); // EXPECTED: [2,4]

	optguardreport ogrep;
	minbleicoptguardresults(state, ogrep);
	//printf("%s\n", ogrep.badgradsuspected ? "true" : "false"); // EXPECTED: false
	//printf("%s\n", ogrep.nonc0suspected ? "true" : "false"); // EXPECTED: false
	//printf("%s\n", ogrep.nonc1suspected ? "true" : "false"); // EXPECTED: false

	sol->updateAt(0, v_init[0]);
	sol->updateAt(1, v_init[1]);

	if (v_init[0] < 10e-9 || v_init[1] < 10e-9) {
		int stop = 1;
		//cerr << "Very small: " << v_init[0] << " / " << v_init[1] << endl;
	}


	/*cout << "OrgFunc: " << org_func(v_init[0], v_init[1], the_info.x_n, the_info.y_n) << endl;

	for (int i = 0; i < the_info.y_n.size(); i++) {
		cout << the_info.x_n[i] << " : " << v_init[0] + v_init[1] * g_x_n(the_info.x_n[i]) << endl;
	}*/

	return 0;
}

int raRSDParameterEstimatorFunctionGradAnalytical(
	const double& v0, const double& v1,
	const RMSEParameterEstimatorAlglibParameter& params,
	double& f, double* grad)
{
	f = 0;
	grad[0] = 0;
	grad[1] = 0;

	for (int i = 0; i < params.y_n.size(); i++)
	{
		double x_i = params.x_n[i];
		double y_i = params.y_n[i];
		double gxn = params.sol->evaluateConstantTermAt(x_i);

		// RMSE formula, but squared
		f += ((y_i - v0 - v1 * gxn) / y_i) * ((y_i - v0 - v1 * gxn) / y_i);

		// RMSE-grad0 formula
		grad[0] += (-2.0 / y_i) * ((y_i - v0 - v1 * gxn) / y_i);
		grad[1] += (-2.0 / y_i) * gxn * ((y_i - v0 - v1 * gxn) / y_i);
	}

	return 0;
}

B<double> RMSEParameterEstimatorFunction(const B<double>& v0, const B<double>& v1,
	const RMSEParameterEstimatorAlglibParameter& params) {
	B<double> z = 0;

	for (int i = 0; i < params.y_n.size(); i++)
	{
		double x_i = params.x_n[i];
		double y_i = params.y_n[i];
		double gxn = params.sol->evaluateConstantTermAt(x_i);

		// RMSE formula
		z += ((y_i - v0 - v1 * gxn)) * ((y_i - v0 - v1 * gxn));		
	}
	z /= (double)params.y_n.size();
	z = sqrt(z);
	return z;
}


void RMSEParameterEstimatorFunctionGradEval(const real_1d_array& v_i, double& func,
	real_1d_array& grad, void* ptr)
{

	RMSEParameterEstimatorAlglibParameter& params =
		*((RMSEParameterEstimatorAlglibParameter*)ptr);

	B<double> v0 = v_i[0];
	B<double> v1 = v_i[1];
	B<double> fval;

	fval = RMSEParameterEstimatorFunction(v0, v1, params);
	fval.diff(0, 1);        // Differentiate f (index 0 of 1)

	func = fval.x();
	grad[0] = v0.d(0); // Value of df/dx (index 0 of 1)
	grad[1] = v1.d(0); // Value of df/dy (index 0 of 1)

	//double f1 = 0;
	//double grad2[2];
}