#include<stdlib.h>
#include<iostream>
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<climits>
#include<list>
#include<algorithm>
#include<numeric>
#include<vector>
#include<iostream>
#include<random>
#include<chrono>
#include "BSM.hpp"
#include "Discretization_schemes.h"

using namespace std;

int main() {

	double S0 = 100;
	double K = 110;
	double r = 0.05;
	double v = 0.2;
	double T = 1;
	unsigned num_intervals = 250;
	int num_sims = 1000;

//	cout << "\nBSM Price Put option = " << BSM_Put(S0, r, v, T, K);
//	cout << "\nBSM Price Call option = " << BSM_Call(S0, r, v, T, K);

	vector<double> euler_spot_prices(num_intervals, S0);
	vector<double> milstein_spot_prices(num_intervals, S0);

	double sum_put_euler = 0;
	double sum_put_milstein = 0;
	double sum_call_euler = 0;
	double sum_call_milstein = 0;

	for (int n = 0; n < num_sims; n++) {

		GBM_EULER(euler_spot_prices, r, v, T);
		GBM_MILSTEIN(milstein_spot_prices, r, v, T);

		sum_put_euler += std::max(K - euler_spot_prices[num_intervals - 1], 0.0);
		sum_put_milstein += std::max(K - milstein_spot_prices[num_intervals - 1], 0.0);	
//		sum_call_euler += std::max(euler_spot_prices[num_intervals - 1] - K, 0.0);
//		sum_call_milstein += std::max(milstein_spot_prices[num_intervals - 1] - K, 0.0);
	}

	double put_euler = sum_put_euler / static_cast<double>(num_sims);
	double put_milstein = sum_put_milstein / static_cast<double>(num_sims);
//	double call_euler = sum_call_euler / static_cast<double>(num_sims);
//	double call_milstein = sum_call_milstein / static_cast<double>(num_sims);


//	cout << "\n\nPrice Put option with Euler = " << put_euler;
//	cout << "\nPrice Put option with Milstein = " << put_milstein;
//	cout << "\n\nPrice Call option with Euler = " << call_euler;
//	cout << "\nPrice Call option with Milstein = " << call_milstein;
//	cout << "\n";

	int L = 10 + 1; // nbr of levels
	int M = 100000; //pow(10.0,5.0);
	vector<double> M_l(L);
	double mu = 0;
	double sigma = 0;
//	vector<double> X_f(S0);
//	vector<double> X_c;
	double M0 = 100;
	vector<double> X(M0,S0);

		for (int l = 0; l < L; l++) {
		M_l[l] = std::pow(2.0,l);
	
		if (l == 0) {
			GBM_EULER(X, r, v, T);
		}
		else {
			for (int n = 0; n < static_cast<int>(M_l[l]); n++) {
				vector<double> X_f(M_l[l],S0);
				vector<double> X_c(M_l[l], S0);
				GBM_EULER_V2(X_f, X_c, M_l[l], r, v, T); // joint
				//mu += std::max(X_f[M_l[l]-1] - X_c[M_l[l] - 1] - K,0.0);
				//sigma += sigma;
			}
			//double mean = mu / static_cast<double>(M);
		}
	}

		cout << M_l[L-1];
//		cout << X_f[L - 1];

	return 0;
}

