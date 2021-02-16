#include "raRSDParameterEstimator.h"
#include "EigenParameterEstimator.h"
#include <iostream>

raRSDParameterEstimator::raRSDParameterEstimator()
{
	_mdb = nullptr;
}

raRSDParameterEstimator::raRSDParameterEstimator(MeasurementDB* mdb) : _mdb(mdb)
{
}

raRSDParameterEstimator::~raRSDParameterEstimator()
{
	_mdb = nullptr;
}

int raRSDParameterEstimator::estimateParameters(AbstractSolution* sol, double newrelerr) {
	int m = _mdb->get_size();
	int n = 2; // Nummber of variables

	raRSDParameterEstimatorAlglibParameter params;

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
	if (eigpar.estimateParameters(sol))
		return ERR_GENERAL_INVALID;

	double par1 = sol->getAt(0);
	double par2 = sol->getAt(1);

	real_1d_array v_init; v_init.setlength(2);
	v_init[0] = par1; v_init[1] = par2;
	real_1d_array s = "[1 , 1]"; // scale does not make difference
	real_2d_array c = "[[1, 0, 0.1], [0, 1 , 10e-9]]";
	integer_1d_array ct = "[1,1]"; // >, >

	minbleicstate state;
	double epsg = 0;
	double epsf = 0;
	double epsx = 0.0000001;
	ae_int_t maxits = 0;

	minbleiccreate(v_init, state);
	minbleicsetlc(state, c, ct);
	minbleicsetscale(state, s);
	minbleicsetcond(state, epsg, epsf, epsx, maxits);

	minbleicoptguardsmoothness(state);
	minbleicoptguardgradient(state, 0.001);
	minbleicreport rep;

	try
	{ 
		alglib::minbleicoptimize(state, raRSDParameterEstimatorFunctionGradEval, nullptr, &params);
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

	if (v_init[0] < 10e-9 || v_init[0] < 10e-9) {
		int stop = 1;
		cin >> stop;
	}
		

	/*cout << "OrgFunc: " << org_func(v_init[0], v_init[1], the_info.x_n, the_info.y_n) << endl;

	for (int i = 0; i < the_info.y_n.size(); i++) {
		cout << the_info.x_n[i] << " : " << v_init[0] + v_init[1] * g_x_n(the_info.x_n[i]) << endl;
	}*/

	return 0;
}

int raRSDParameterEstimatorFunctionGradAnalytical(
	const double& v0, const double& v1,
	const raRSDParameterEstimatorAlglibParameter& params,
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

		// nnrRSD formula, but squared
		f += ((y_i - v0 - v1 * gxn) / y_i) * ((y_i - v0 - v1 * gxn) / y_i);

		// nnrRSD-grad0 formula
		grad[0] += (-2.0 / y_i) * ((y_i - v0 - v1 * gxn) / y_i);
		grad[1] += (-2.0 / y_i) * gxn * ((y_i - v0 - v1 * gxn) / y_i);
	}

	return 0;
}

B<double> raRSDParameterEstimatorFunction(const B<double>& v0, const B<double>& v1,
	const raRSDParameterEstimatorAlglibParameter& params) {
	B<double> z = 0;

	for (int i = 0; i < params.y_n.size(); i++)
	{
		double x_i = params.x_n[i];
		double y_i = params.y_n[i];
		double gxn = params.sol->evaluateConstantTermAt(x_i);

		// nnrRSD formula, but squared
		z += ((y_i - v0 - v1 * gxn) / y_i) * ((y_i - v0 - v1 * gxn) / y_i);
	}

	return z;
}


void raRSDParameterEstimatorFunctionGradEval(const real_1d_array& v_i, double& func,
	real_1d_array& grad, void* ptr)
{
	//
	// this callback calculates f(x0,x1) = 100*(x0+3)^4 + (x1-3)^4
	// and its derivatives df/d0 and df/dx1
	//

	raRSDParameterEstimatorAlglibParameter& params =
		*((raRSDParameterEstimatorAlglibParameter*)ptr);

	B<double> v0 = v_i[0];
	B<double> v1 = v_i[1];
	B<double> fval;

	/*fval = raRSDParameterEstimatorFunction(v0, v1, params);
	fval.diff(0, 1);        // Differentiate f (index 0 of 1)

	func = fval.x();
	grad[0] = v0.d(0); // Value of df/dx (index 0 of 1)
	grad[1] = v1.d(0); // Value of df/dy (index 0 of 1)


	double f1 = 0;
	double grad2[2];*/

	raRSDParameterEstimatorFunctionGradAnalytical(v0.x(), v1.x(), params, func, &grad[0]);
}