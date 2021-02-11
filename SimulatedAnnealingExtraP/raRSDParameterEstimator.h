#ifndef ANRSDESTIMATOR_H
#define ANRSDESTIMATOR_H

#include "MeasurementDB.h"
#include "ParameterEstimatorInterface.h"
#include "AbstractSolution.h"
#include "alglib/optimization.h"
#include "FADBAD++/badiff.h"
#include <vector>

using namespace alglib;
using namespace std;
using namespace fadbad;

typedef struct  {
	MeasurementDB* mdb;
	AbstractSolution* sol;
	vector<double> x_n;
	vector<double> y_n;
} raRSDParameterEstimatorAlglibParameter;

void raRSDParameterEstimatorFunctionGradEval(const real_1d_array& v_i, double& func,
	real_1d_array& grad, void* ptr);

B<double> raRSDParameterEstimatorFunction(const B<double>& v0, const B<double>& v1,
	const raRSDParameterEstimatorAlglibParameter& params);

class raRSDParameterEstimator : public ParameterEstimatorInterface {
public:
	raRSDParameterEstimator();
	raRSDParameterEstimator(MeasurementDB* mdb);
	~raRSDParameterEstimator();
	void estimateParameters(AbstractSolution* sol, double newrelerr=0.0);

	int raRSDParameterEstimatorFunctionGradAnalytical(const double& v0, const double& v1,
		const raRSDParameterEstimatorAlglibParameter& params, double& f, double* grad);

private:
	MeasurementDB* _mdb;
};

#endif // ! NNRRSDESTIMATOR_H

