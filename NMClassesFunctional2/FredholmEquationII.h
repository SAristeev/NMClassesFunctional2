#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <functional>
#include <iomanip>
#include "functions.h"

// The class contains the task of completely
// All functions are stored inside the class

class FredholmEquationII {
public:
	int nodes; // count of nodes of unit integration segment
	int p; // order of accuracy
	int var; // test case number
	double a, b; // equation boundaries
	double h; // uniform gird size
	double I = 0; // value of the solution integral. It is initialized only after solving the equation

	std::function<double(double)> f; // right - hand side of the equation
	std::function<double(double, double)> K; // kernel
	std::function<double(double)> u = std::function<double(double)>([](double x) {return 0; }); // solution functor. Zero initialized by default
	// u(x) = sum (G(i)*K(x,i)*U[i] + f(x))
	std::vector<double> U;
	std::vector<double> G; // coefficients to calculate the integral on integration segment
	std::set<double> mesh; // uneven grid

	FredholmEquationII(double a, double b, std::function<double(double)> f, std::function<double(double, double)> K, int nodes, int var);
	void solve();

	double computeIntegral(); // compute integral of solved function u. Used adaptive mesh

	void PlotAdaptiveMesh(); // can be used after computeIntegral
	void Plot(int n); // plot by n points

};

