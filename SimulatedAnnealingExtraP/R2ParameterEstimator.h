#ifndef R2PARAMETERESTIMATOR_H
#define R2PARAMETERESTIMATOR_H

#include "MeasurementDB.h"
#include "ParameterEstimatorInterface.h"
#include "AbstractSolution.h"
#include "alglib/optimization.h"
#include "FADBAD++/badiff.h"
#include <vector>

using namespace alglib;
using namespace std;
using namespace fadbad;

typedef struct {
	MeasurementDB* mdb;
	AbstractSolution* sol;
	vector<double> x_n;
	vector<double> y_n;
} R2ParameterEstimatorAlglibParameter;

void R2ParameterEstimatorFunctionGradEval(const real_1d_array& v_i, double& func,
	real_1d_array& grad, void* ptr);

B<double> R2ParameterEstimatorFunction(const B<double>& v0, const B<double>& v1,
	const R2ParameterEstimatorAlglibParameter& params);

class R2ParameterEstimator : public ParameterEstimatorInterface {
public:
	R2ParameterEstimator();
	R2ParameterEstimator(MeasurementDB* mdb);
	~R2ParameterEstimator();
	int estimateParameters(AbstractSolution* sol, double newrelerr = 0.0);

	int R2ParameterEstimatorFunctionGradAnalytical(const double& v0, const double& v1,
		const R2ParameterEstimatorAlglibParameter& params, double& f, double* grad);

private:
	MeasurementDB* _mdb;
};

#endif
