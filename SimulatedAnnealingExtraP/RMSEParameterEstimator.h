#ifndef RMSERESTIMATOR_H
#define RMSERESTIMATOR_H

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
} RMSEParameterEstimatorAlglibParameter;

void RMSEParameterEstimatorFunctionGradEval(const real_1d_array& v_i, double& func,
	real_1d_array& grad, void* ptr);

B<double> RMSEParameterEstimatorFunction(const B<double>& v0, const B<double>& v1,
	const RMSEParameterEstimatorAlglibParameter& params);

class RMSEParameterEstimator : public ParameterEstimatorInterface {
public:
	RMSEParameterEstimator();
	RMSEParameterEstimator(MeasurementDB* mdb);
	~RMSEParameterEstimator();
	int estimateParameters(AbstractSolution* sol, double newrelerr = 0.0);

	int RMSEParameterEstimatorFunctionGradAnalytical(const double& v0, const double& v1,
		const RMSEParameterEstimatorAlglibParameter& params, double& f, double* grad);

private:
	MeasurementDB* _mdb;
};

#endif // ! NNRRSDESTIMATOR_H

