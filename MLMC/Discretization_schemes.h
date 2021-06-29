#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>

using namespace std;

double gaussian_box_muller() {
    double x = 0.0; double y = 0.0; double euclid_sq = 0.0;

    do {
        x = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
        y = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
        euclid_sq = x * x + y * y;
    } while (euclid_sq >= 1.0);

    return x * sqrt(-2 * log(euclid_sq) / euclid_sq);
}


void GBM_EULER(std::vector<double>& spot_prices, const double& r, const double& v, const double& T) { 

    double dt = T / static_cast<double>(spot_prices.size());
    double drift = exp(dt * (r - 0.5 * v * v));
    double vol = sqrt(v * v * dt);

    for (int i = 1; i < spot_prices.size(); i++) {
        double gauss_bm = gaussian_box_muller();
        spot_prices[i] = spot_prices[i - 1] * drift * exp(vol * gauss_bm);
    }
}

void GBM_EULER_V2(std::vector<double>& X_f, std::vector<double>& X_c, const double& time_steps, const double& r, const double& v, const double& T) {

    double dt = T / time_steps;

    for (int i = 1; i < time_steps; i++) {
        double gauss_bm1 = gaussian_box_muller();
        X_f[i] = X_f[i - 1] + r * dt + v * (sqrt(dt) * gauss_bm1);
        double gauss_bm2 = gaussian_box_muller();
       // X_c[i] = X_c[i - 2] + r * (T / (static_cast<double>(time_steps) - 1)) + v * (sqrt(dt) * (gauss_bm1 + gauss_bm2));
    }
}

void GBM_MILSTEIN(std::vector<double>& spot_prices, const double& r, const double& v, const double& T) {

    double dt = T / static_cast<double>(spot_prices.size());

    for (int t = 1; t < spot_prices.size(); t++) {
        double gauss_bm = gaussian_box_muller();
        spot_prices[t] = spot_prices[t - 1] + r * dt + v * (sqrt(dt)*gauss_bm) + 0.5 * v * v * ((sqrt(dt) * gauss_bm)* (sqrt(dt) * gauss_bm) - dt);
    }
}
