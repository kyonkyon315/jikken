#include "solve2.h"
#include <iostream>
#include <string>
//[cgs]
#define pi 3.14159265359
#define a0 0.529e-8
#define e  4.803204673e-10
#define k 1.381e-16
#define me 9.1093837e-28
#define h 6.626e-27
#define eV 1.6022e-12//[erg]
#define sigma 5.6704e-5
#define aumstrong 1e-8;

double HiEnergy(size_t n) {
	if (n == 0)throw 0;
	return -pow(e, 2.) / (2. * a0 * pow((double)n, 2.));
}

double BoltzmannFactor(double T,size_t n) {
	return 2. * pow((double)n, 2.) * exp(-HiEnergy(n) / (k * T));
}

double Z(double T) {
	double ans = 0.;
	for (int i = 1; i < 11; ++i) {
		ans += BoltzmannFactor(T,i);
		//2 n^2 exp{- Ei/kT}
	}
	return ans;
}

double f(double T) {
	const double kai = 13.6 * eV;
	double ans = 1. / (Z(T) / 2.);
	ans *= pow(2*pi*me*k*T/pow(h,2.), 1.5);
	ans *= exp(-kai / (k * T));
	if (ans < 0.)std::cout << "error" << std::flush;
	return ans;
}

double A(double Nh, double T) {
	return Nh * f(T);
}

double x_ratio(double Nh, double T) {
	double A_value = A(Nh, T);
	if (A_value > 1e5)return 1.-2/A_value;
	return (-A_value + sqrt(pow(A_value,2.) + 4. * A_value)) / 2.;
	//return A_value * (-1. + sqrt(1. + 4. / A_value)) / 2.;
}

double h_func(double T) {
	const double kai = 0.75 * eV;
	double ans = 4.;
	ans *= pow(2 * pi * me * k * T / pow(h, 2.), 1.5);
	ans *= exp(-kai / (k * T));
	return ans;
}

double Nminus_ratio(double Nh, double T) {
	double Nplus_ratio = x_ratio(Nh, T);
	return (1 - Nplus_ratio) * Nplus_ratio * Nh / h_func(T);
}

double xVSa(double a) {
	double Nh = 2.4e-7 / pow(a, 3.);
	double T = 2.7 / a;

	double A_value = A(Nh, T);

	if (A_value > 1e10)return 1 - 2 / A_value;
	return (-A_value + sqrt(pow(A_value, 2.) + 4. * A_value)) / 2.;
}

double sigmaBFn(double lambda, size_t n) {
	const double alpha0 = 1.0449e-26;
	const double g_bf = 1.;
	const double R = 1.0968e-3;

	if (lambda < pow((double)n, 2.) / R) {
		return alpha0 * g_bf * pow(lambda, 3.) / pow((double)n, 5.);
	}
	else {
		return 0.;
	}
}

double sigmaFFminus(double lambda) {
	const double a0ff = 1.99654;
	const double a1ff = -1.18267e-5;
	const double a2ff = 2.64243e-6;
	const double a3ff = -4.40524e-10;
	const double a4ff = 3.23992e-14;
	const double a5ff = -1.39568e-18;
	const double a6ff = 2.78701e-23;

	return (a0ff + a1ff * lambda
		+ a2ff * pow(lambda, 2.) + a3ff * pow(lambda, 3.)
		+ a4ff * pow(lambda, 4.) + a5ff * pow(lambda, 5.)
		+ a6ff * pow(lambda, 6.)) * 1e-18;
}
double sigmaBFminus(double lambda,double T,double x) {
	double theta = 5040. / T;
	const double a00bf = -2.2763, a01bf = -1.6850, a02bf = 0.76661, a03bf = -0.053346;
	const double a10bf = 15.2827, a11bf = -9.2846, a12bf = 1.99381, a13bf = -0.142631;
	const double a20bf = -197.789, a21bf = 190.266, a22bf = -67.9775, a23bf = 10.6913,a24bf=-0.625151;

	double f0 = a00bf + a01bf * log10(lambda) + a02bf * pow(log10(lambda), 2.) + a03bf * pow(log10(lambda), 3.);
	double f1 = a10bf + a11bf * log10(lambda) + a12bf * pow(log10(lambda), 2.) + a13bf * pow(log10(lambda), 3.);
	double f2 = a20bf + a21bf * log10(lambda) + a22bf * pow(log10(lambda), 2.) + a23bf * pow(log10(lambda), 3.)
		       +a24bf*pow(log10(lambda),4.);

	double Pe = x * k * T;

	return (Pe * pow(10., f0 + f1 * log(theta) + f2 * pow(log(theta), 2.))) * 1e-26;
}

void alpha_lambda(double T, double Nh,size_t label) {
	std::vector<double> lambdas;
	for (double lambda = 3000.; lambda <= 20000.; lambda += 100.) {
		lambdas.push_back(lambda);
	}

	//N(H(N))/N(H)
	std::vector<double> N_HnvsN_H(10);
	double Z_ = Z(T);
	for (size_t i = 0; i < N_HnvsN_H.size(); i++) {
		N_HnvsN_H[i] = BoltzmannFactor(T, i + 1) / Z_;
	}
	
	double x = x_ratio(Nh, T);
	//N(H)/N_H
	double N_HvsNh = 1 - x;
	
	//N(H^-)/N_H
	double NminusvsNh=Nminus_ratio(Nh,T);

	std::vector<double> alpha_lambdas(lambdas.size());
	for (size_t i = 0; i < lambdas.size(); i++) {
		double sum = 0.;
		for (size_t j = 0; j < 10; j++) {
			sum += sigmaBFn(lambdas[i], j + 1) * N_HnvsN_H[j];
		}
		alpha_lambdas[i] = (N_HvsNh *( sum + sigmaFFminus(lambdas[i])) + NminusvsNh * sigmaBFminus(lambdas[i],T,x)) * Nh;
	}

	savetxt(lambdas, alpha_lambdas,"question5_fig11" + std::to_string(label) + ".txt");
}

void solve2() {

	std::vector<double>Ts;
	std::vector<double>Zs;
	for (double T = 1000; T <= 100000; T += 100) {
		Ts.push_back(T);
		Zs.push_back(Z(T));
	}
	savetxt(Ts, Zs, "question2_fig1.txt");

	std::vector<std::vector<double>>Nratio(10,std::vector<double>(Ts.size()));
	for (size_t i = 0; i < Nratio.size(); i++) {
		for (size_t j = 0; j < Ts.size(); j++) {
			Nratio[i][j] = BoltzmannFactor(Ts[j], i + 1) / Zs[j];
		}
		savetxt(Ts, Nratio[i], "question2_fig2"+std::to_string(i+1)+".txt");
	}

	std::vector<double>Nhs = { 1,1e10,1e20 };
	std::vector<std::vector<double>>x_ratioVS_T(Nhs.size(), std::vector<double>(Ts.size()));

	for (size_t i = 0; i < Nhs.size(); i++) {
		for (size_t j = 0; j < Ts.size(); j++) {
			x_ratioVS_T[i][j] = x_ratio(Nhs[i], Ts[j]);
		}
		savetxt(Ts, x_ratioVS_T[i], "question2_fig3" + std::to_string(i) + ".txt");
	}

	std::vector<double>Ts_minus;
	std::vector<std::vector<double>> Nminus_ratioVS_T(Nhs.size());

	for (double T = 1000; T <= 10000; T += 100) {
		Ts_minus.push_back(T);
		for (size_t i = 0; i < Nhs.size(); i++) {
			Nminus_ratioVS_T[i].push_back(Nminus_ratio(Nhs[i], T));
		}
	}
	for (size_t i = 0; i < Nhs.size(); i++) {
		savetxt(Ts_minus, Nminus_ratioVS_T[i], "question2_fig4" + std::to_string(i) + ".txt");
	}


	std::vector<double > as;
	std::vector<double>xVSas;
	for (double a = 1e-4; a <= 1e-2; a += 1e-4) {
		as.push_back(a);
		xVSas.push_back(xVSa(a));
	}
	savetxt(as, xVSas, "question2_fig5.txt");

	std::vector<std::pair<double, double>>inputs = {
		{5000,1e18},
		{6000,1e16},
		{8000,1e14},
		{12000,1e13},
	};

	for(size_t i=0;i<inputs.size();i++){
		alpha_lambda(inputs[i].first, inputs[i].second, i + 1);
	}
}

