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
#include "Payoffs.h"

using namespace std;

int main() {

	double S0 = 90;
	double K = 110;
	double r = 0.05;
	double v = 0.2;
	double T = 1;
	double act = exp(-r * T);
	double B = 100;
	unsigned num_intervals = 250;
	int num_sims = 10000;

	cout << "\nFor S0 = " << S0 << ", K = " << K << ", r = " << r << ", v = " << v << " and B = " << B << " : ";

	cout << "\n\nBSM Put price = " << BSM_Put(S0, r, v, T, K);
	cout << "\nBSM Call price = " << BSM_Call(S0, r, v, T, K);

	cout << "\n\nBSM Binary Put price = " << BSM_BinPut(S0, r, v, T, K);
	cout << "\nBSM Binary Call price = " << BSM_BinCall(S0, r, v, T, K);

	cout << "\n\nBSM Up-and-Out Put price = " << BSM_UOP(S0, r, v, T, K, B);
	cout << "\nBSM Up-and-Out Call price = " << BSM_UOC(S0, r, v, T, K, B);

	cout << "\n\nBSM Up-and-In Put price = " << BSM_UIP(S0, r, v, T, K, B);
	cout << "\nBSM Up-and-In Call price = " << BSM_UIC(S0, r, v, T, K, B);

	cout << "\n\nBSM Down-and-Out Put price = " << BSM_DOP(S0, r, v, T, K, B);
	cout << "\nBSM Down-and-Out Call price = " << BSM_DOC(S0, r, v, T, K, B);

	cout << "\n\nBSM Down-and-In Put price = " << BSM_DIP(S0, r, v, T, K, B);
	cout << "\nBSM Down-and-In Call price = " << BSM_DIC(S0, r, v, T, K, B);
	
	vector<double> euler_spot_prices(num_intervals, S0);
	vector<double> milstein_spot_prices(num_intervals, S0);

	double put_euler = 0; double call_euler = 0;
	double put_milstein = 0; double call_milstein = 0;
	double put_euler_squared = 0; double call_euler_squared = 0;
	double proba_crossing = 0;
	double BB_MC_UOP = 0; double BB_MC_UIP = 0;	double BB_MC_DOP = 0; double BB_MC_DIP = 0;
	double BB_MC_UOC = 0; double BB_MC_UIC = 0;	double BB_MC_DOC = 0; double BB_MC_DIC = 0;
	double MC_UIC = 0; double MC_UIP = 0; double MC_UOC = 0; double MC_UOP = 0;
	double MC_DIC = 0; double MC_DIP = 0; double MC_DOC = 0; double MC_DOP = 0;

	double put_heston = 0; double call_heston = 0; ///
	double V0 = 0.225;
	vector<double> vol_paths(num_intervals, V0);
	double kappa = 1;
	double theta = 0.25;
	double xi = 0.02;
	double rho = 0.5;

	for (int n = 0; n < num_sims; n++) {

		GBM_EULER(euler_spot_prices, r, v, T);
		//GBM_MILSTEIN(milstein_spot_prices, r, v, T);

		put_euler += act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
		call_euler += act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
		put_euler_squared += (pow(Put(euler_spot_prices[num_intervals - 1], K), 2)) / (static_cast<double>(num_sims) - 1);
		call_euler_squared += (pow(Call(euler_spot_prices[num_intervals - 1], K), 2)) / (static_cast<double>(num_sims) - 1);

/*		HESTON_EULER(euler_spot_prices, vol_paths, r, v, T, kappa, theta, xi, rho);
		put_heston += act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);//
		call_heston += act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);//
*/
		if (*max_element(euler_spot_prices.begin(), euler_spot_prices.end()) >= B) {
			MC_UIC += act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_UIP += act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_UOC += 0; MC_UOP += 0;
		}
		else {
			MC_UIC += 0; MC_UIP += 0;
			MC_UOC += act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_UOP += act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
		}

		if (*min_element(euler_spot_prices.begin(), euler_spot_prices.end()) <= B) {
			MC_DIC += act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_DIP += act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);;
			MC_DOC += 0; MC_DOP += 0;
		}
		else {
			MC_DOC += act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_DOP += act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_DIC += 0; MC_DIP += 0;
		}

		vector<double> proba_up(num_intervals - 1, 1); vector<double> proba_down(num_intervals - 1, 1);
		double proba_up_cross = 1; double proba_down_cross = 1;

		Brownian_Bridge(euler_spot_prices, proba_up, proba_down, r, v, T, B, proba_up_cross, proba_down_cross);

		// proba_crossing = (std::accumulate(proba.begin(), proba.end(), proba[0], std::multiplies<double>()));

		BB_MC_UOP += act * Put(euler_spot_prices[num_intervals - 1], K) * (proba_up_cross) / static_cast<double>(num_sims);
		BB_MC_UOC += act * Call(euler_spot_prices[num_intervals - 1], K) * (proba_up_cross) / static_cast<double>(num_sims);

		BB_MC_UIP += act * Put(euler_spot_prices[num_intervals - 1], K) * (1 - proba_up_cross) / static_cast<double>(num_sims);
		BB_MC_UIC += act * Call(euler_spot_prices[num_intervals - 1], K) * (1 - proba_up_cross) / static_cast<double>(num_sims);

		BB_MC_DOP += act * Put(euler_spot_prices[num_intervals - 1], K) * (proba_down_cross) / static_cast<double>(num_sims);
		BB_MC_DOC += act * Call(euler_spot_prices[num_intervals - 1], K) * (proba_down_cross) / static_cast<double>(num_sims);

		BB_MC_DIP += act * Put(euler_spot_prices[num_intervals - 1], K) * (1 - proba_down_cross) / static_cast<double>(num_sims);
		BB_MC_DIC += act * Call(euler_spot_prices[num_intervals - 1], K) * (1 - proba_down_cross) / static_cast<double>(num_sims);

	}
	
	double put_MC_std = pow((put_euler_squared - pow(put_euler, 2)), 0.5);
	double call_MC_std = pow((call_euler_squared - pow(call_euler, 2)), 0.5);

	cout << "\n\n\nStandard MC with nbr of simulations = " << num_sims << " :";

	cout << "\n\nMC Put price GBM Euler = " << put_euler << " with MC std dev = " << put_MC_std * 100 << "%";
	cout << " and MC std error = " << put_MC_std / (static_cast<double>(num_sims)) * 100 << "%";

	cout << "\nMC Call price Euler = " << call_euler << " with MC std dev = " << call_MC_std * 100 << "%";
	cout << " and MC std error = " << call_MC_std / (static_cast<double>(num_sims)) * 100 << "%";

/*	cout << "\n\nMC Put price Heston Euler = " << put_heston;
	cout << "\nMC Call price Heston Euler = " << call_heston;
*/
	cout << "\n\n\nMC Up-and-Out Put price Euler = " << MC_UOP;
	cout << "\nMC Up-and-Out Call price Euler = " << MC_UOC;

	cout << "\n\nMC Up-and-In Put price Euler = " << MC_UIP;
	cout << "\nMC Up-and-In Call price Euler = " << MC_UIC;

	cout << "\n\nMC Down-and-Out Put price Euler = " << MC_DOP;
	cout << "\nMC Down-and-Out Call price Euler = " << MC_DOC;

	cout << "\n\nMC Down-and-In Put price Euler = " << MC_DIP;
	cout << "\nMC Down-and-In Call price Euler = " << MC_DIC;

// Brownian Bridge

	cout << "\n\n\nBrownian Bridge MC :";
	cout << "\n\nMC BB Up-and-Out Put price Euler = " << BB_MC_UOP;
	cout << "\nMC BB Up-and-Out Call price Euler = " << BB_MC_UOC;

	cout << "\n\nMC BB Up-and-In Put price Euler = " << BB_MC_UIP;
	cout << "\nMC BB Up-and-In Call price Euler = " << BB_MC_UIC;

	cout << "\n\nMC BB Down-and-Out Put price Euler = " << BB_MC_DOP;
	cout << "\nMC BB Down-and-Out Call price Euler = " << BB_MC_DOC;

	cout << "\n\nMC BB Down-and-In Put price Euler = " << BB_MC_DIP;
	cout << "\nMC BB Down-and-In Call price Euler = " << BB_MC_DIC;

//	cout << "\nMC Put price with Milstein = " << put_milstein;
//	cout << "\nPrice Call option with Milstein = " << call_milstein;
//	cout << "\n";

// MLMC

	int L = 10 + 1; // nbr of levels
//	int M = 100000; //pow(10.0,5.0);
	vector<double> N_l(L);
	double mu = 0;
	double mu0 = 0;
	double sigma = 0;
	double sigma0 = 0;
	double N0 = 1000; // std::pow(2.0, L + 1);
	vector<double> X(N0,S0);

		for (int l = 0; l < L; l++) {
		N_l[l] = std::pow(2.0,l);
	
		if (l == 0) {
			for (int n = 0; n < N0; n++) {
				GBM_EULER(X, r, v, T);
				mu0 += act * Put(X[N0 - 1], K) / N0;
				sigma0 += (pow(Put(X[N0 - 1], K),2))/ N0;
			}
		}
		else {
			for (int n = 0; n < static_cast<int>(N_l[l]); n++) {

				vector<double> X_f(N_l[l],S0);
				vector<double> X_c(N_l[l], S0);

				GBM_EULER_V2(X_f, X_c, N_l[l], r, v, T); ///

				mu += exp(-r * T / N_l[l]) * (Put(X_f[N_l[l] - 1], K) - Put(X_c[N_l[l - 1] - 1], K)) / N_l[l]; // different time steps
				sigma += (pow(Put(X_f[N_l[l] - 1], K) - Put(X_c[N_l[l - 1] - 1], K), 2)) / (N_l[l] - 1);
			}
		}
	}

		double put_MLMC_std = pow(sigma - pow(mu, 2) + sigma0 - pow(mu0, 2), 0.5);

		cout << "\n\n\nMLMC :";
		cout << "\n\nMLMC Put price Euler = " << mu + mu0 << " with MLMC std dev = " << put_MLMC_std * 100 << "%";
		cout << " and MLMC std error = " << put_MLMC_std / (std::pow(2.0,L)) * 100 << "%";
		
		cout << "\n\nmu0 = " << mu0;
		cout << "\nmu = " << mu;

		cout << "\n\nMLMC UOC = ";


		cout << "\n";

	return 0;
}

