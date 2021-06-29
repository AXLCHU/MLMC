#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <boost/math/distributions/normal.hpp>

using boost::math::normal;

double BSM_Call(double& S0, const double& r, const double& v, const double& T, const double& K) {

	double d1 = (1 / v * sqrt(T)) * (log(S0 / K) + (r + v * v / 2) * T);
	double d2 = d1 - v * sqrt(T);
	normal s;
	double CallBS = S0 * cdf(s,d1) - K * cdf(s, d2) * exp(-r*T);
	return CallBS;
}

double BSM_Put(double& S0, const double& r, const double& v, const double& T, const double& K) {

	double d1 = (1 / v * sqrt(T)) * (log(S0 / K) + (r - v * v / 2) * T);
	double d2 = d1 - v * sqrt(T);
	normal s;
	double PutBS = K * cdf(s, -d2) * exp(-r * T) - S0 * cdf(s, -d1);
	return PutBS;
}