#include "solve1.h"
#include <iostream>
#include <vector>

static double H(double x) {
	return 1. - exp(-x) - x / 3.;
}

static double solveEq() {
	double L = 0.1;
	double R = 3.;
	double delta = 0.00000001;
	double newValue;
	while (R - L > delta) {
		newValue = (R + L) / 2.;
		if (H(newValue) > 0.) {
			L = newValue;
		}
		else {
			R = newValue;
		}
	}
	return newValue;
}

void solve1() {
	double ans = solveEq();
	std::cout << ans << "\n";
	std::cout << H(ans) << "\n";

	const double pi = 3.1415926535;
	const double c = 299792458.;
	const double h = 6.62607015e-34;
	const double k = 1.380649e-23;

	std::cout << k * ans / h << "\n";

	double F = pow(pi, 5.) * 2. * pow(k, 4.) / (15. * pow(c, 2.) * pow(h, 3.));
	std::cout << F << "\n";
}

