#include "solve4.h"
#include "Printer.h"
#include <utility>
#include <iostream>
#define sigma  5.670374419e-8
double function4(double x) {
	return exp(-x) / x;
}

std::pair<std::vector<double>,std::vector<double>>numericalInt() {

	double plotInterval = 0.01;

	double deltaT = 0.0001;

	double infinityT = 400.;

	double startX = 0.01;
	double endX = 10.;

	size_t numOfData = (size_t)((endX - startX) / plotInterval) + 1ULL;

	std::vector<double>x(numOfData), y(numOfData);

	double ans = 0.;

	double nextX = endX;

	double T = infinityT;
	size_t id = 0;

	for (;;) {
		double simpson = (function4(T)
			+ 4. * function4(T + deltaT / 2.)
			+ function4(T + deltaT));

		ans += simpson;

		if (T <= nextX) {
			y[id] = ans * (deltaT/6.);
			x[id] = T;
			id++;
			nextX -= plotInterval;

			if (nextX < startX)break;
		}
		T -= deltaT;
	}
	x.resize(id);
	y.resize(id);
	return std::pair<std::vector<double>, std::vector<double>>(std::move(x), std::move(y));
}
const double c0 = -0.57721566;
const double c1 = 0.99999193;
const double c2 = -0.24991055;
const double c3 = 0.05519968;
const double c4 = -0.00976004;
const double c5 = 0.00107857;

const double a3 = 8.5733287401;
const double a2 = 18.0590169730;
const double a1 = 8.6347608925;
const double a0 = 0.2677737343;


const double b3 = 9.5733223454;
const double b2 = 25.6329561486;
const double b1 = 21.0996530827;
const double b0 = 3.9584969228;

double approx(double x) {
	if (x <= 1.) {
		return -log(x)
			+ c0
			+ c1 * x
			+ c2 * pow(x, 2.)
			+ c3 * pow(x, 3.)
			+ c4 * pow(x, 4.)
			+ c5 * pow(x, 5.);
	}
	else {
		return (pow(x, 4.) + a3 * pow(x, 3.) + a2 * pow(x, 2.) + a1 * x + a0)*exp(-x) /
			((pow(x, 4.) + b3 * pow(x, 3.) + b2 * pow(x, 2.) + b1 * x + b0)*x);
	}
}

double calcT(double kapper) {
	if (kapper >= 0.1) {
		return 6500. + 1500. * log10(kapper);
	}
	else {
		return 5200. + 200. * log10(kapper);
	}
}

double integratedFunc(double kapper,double kapperDash) {
	return sigma * pow(calcT(kapperDash), 4.) * approx(abs(kapper - kapperDash));
}
double simpson(double a, double b, size_t n,double kapper) {

	double deltaKapperDash = (b-a)/n;
	double ans = 0.;
	for (double KapperDash = a; KapperDash <= b; KapperDash += deltaKapperDash) {
		ans += (integratedFunc(kapper, KapperDash)
			+ 4. * integratedFunc(kapper, KapperDash + deltaKapperDash / 2.)
			+ integratedFunc(kapper, KapperDash + deltaKapperDash));
		
	}
	return ans * deltaKapperDash / 6.;
}
double calc_J(double kapper) {
	if (kapper < 0.01)throw 1;
	double epsilon = 0.001;
	double deltaKapper = 0.001;
	double kapperInfinity = 400.;
	double ans = simpson(0.001, kapper - epsilon, (size_t)((kapper - epsilon-0.001) / deltaKapper), kapper)
		       + simpson(kapper+epsilon,kapperInfinity , (size_t)((kapperInfinity -kapper- epsilon) / deltaKapper), kapper);
	
	ans +=  sigma * pow(calcT(kapper), 4.) * epsilon * (1 - log(epsilon))/log(10.);
	//std::cout << sigma * pow(calcT(kapper), 4.) * epsilon * (1 - log(epsilon)) << " ";

	return ans;
}


void solve4() {
	auto data= numericalInt();
	std::vector<double>x = data.first;
	std::vector<double>calculatedE1 = data.second;
	savetxt(x, calculatedE1, "question4_fig6.txt");

	std::vector<double>approxValue(x.size());
	std::vector<double>approxDiff(x.size());
	std::vector<double>approxDiffRatio(x.size());

	for (int i = 0; i < x.size(); i++) {
		approxValue[i] = approx(x[i]);
		approxDiff[i] = approxValue[i]-calculatedE1[i];
		approxDiffRatio[i] = approxValue[i] / calculatedE1[i]-1.;
	}

	double deltaKapper = 0.1;
	int N = 55;
	std::vector<double>kappers(N);
	std::vector<double>J(N);
	double kapper = 0.01;
	for (int i = 0; i < N; i++) {
		kappers[i] = kapper;
		J[i] = calc_J(kapper);
		kapper *= 1.2;
	}

	savetxt(x, approxValue, "question4_fig7.txt");
	savetxt(x, approxDiff, "question4_fig8.txt");
	savetxt(x, approxDiffRatio, "question4_fig9.txt");
	savetxt(kappers, J, "question4_fig10.txt");
}


