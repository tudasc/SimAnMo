#include "R2ParameterEstimator.h"
#include "EigenParameterEstimator.h"
#include <iostream>

R2ParameterEstimator::R2ParameterEstimator()
{
	_mdb = nullptr;
}

R2ParameterEstimator::R2ParameterEstimator(MeasurementDB* mdb) : _mdb(mdb)
{
	cout << "Created R2ParameterEstimator" << endl;
}

R2ParameterEstimator::~R2ParameterEstimator()
{
	_mdb = nullptr;
}

int R2ParameterEstimator::estimateParameters(AbstractSolution* sol, double newrelerr) {
	int m = _mdb->get_size();
	int n = 2; // Nummber of variables

	R2ParameterEstimatorAlglibParameter params;

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
		alglib::minbleicoptimize(state, R2ParameterEstimatorFunctionGradEval, nullptr, &params);
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

	return 0;
}

int R2ParameterEstimatorFunctionGradAnalytical(
	const double& v0, const double& v1,
	const R2ParameterEstimatorAlglibParameter& params,
	double& f, double* grad)
{
	cerr << "Calling non-implemted function R2ParameterEstimatorFunctionGradAnalytical" << endl;
	exit(11);

	/*f = 0;
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

	return 0;*/

}

B<double> R2ParameterEstimatorFunction(const B<double>& v0, const B<double>& v1,
	const R2ParameterEstimatorAlglibParameter& params) {
	B<double> r2 = 0;

	double xtimesy = 0.0;
	double xx = 0.0;
	double yy = 0.0;
	double sx = 0.0;
	double sy = 0.0;

	const double n = (double)params.x_n.size();


	for (int i = 0; i < params.y_n.size(); i++)
	{
		std::pair<double, double> act_point{ params.x_n[i] , params.y_n[i] };
		double gxn = params.sol->evaluateConstantTermAt(act_point.first);

		// R2 intermediate data
		sx += act_point.second;
		sy += gxn;
		xx += act_point.second * act_point.second;
		yy += gxn * gxn;
		xtimesy += act_point.second * gxn;
	}

	double nSxy = n * xtimesy;
	double nominator = nSxy - sx * sy;
	double denominator = sqrt((n * xx - sx * sx) * (n * yy - sy * sy));
	double r = nominator / denominator;

	r2 = 1.0 - r *r ; // In order to maximize, not to minimize
	return r2;
}


void R2ParameterEstimatorFunctionGradEval(const real_1d_array& v_i, double& func,
	real_1d_array& grad, void* ptr)
{

	R2ParameterEstimatorAlglibParameter& params =
		*((R2ParameterEstimatorAlglibParameter*)ptr);

	B<double> v0 = v_i[0];
	B<double> v1 = v_i[1];
	B<double> fval;

	fval = R2ParameterEstimatorFunction(v0, v1, params);
	fval.diff(0, 1);        // Differentiate f (index 0 of 1)

	func = fval.x();
	grad[0] = v0.d(0); // Value of df/dx (index 0 of 1)
	grad[1] = v1.d(0); // Value of df/dy (index 0 of 1)

	//double f1 = 0;
	//double grad2[2];
}