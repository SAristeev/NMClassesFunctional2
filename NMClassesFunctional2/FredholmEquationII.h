#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <functional>
#include <iomanip>
#include "functions.h"

class FredholmEquationII {
public:
	int nodes, p, var;
	double a, b, h, I;
	std::function<double(double)> f;
	std::function<double(double, double)> K;
	std::function<double(double)> u = std::function<double(double)>([](double x) {return 0; });
	std::vector<double> U;
	std::vector<double> G;
	std::set<double> mesh;
	FredholmEquationII(double a, double b, std::function<double(double)> f, std::function<double(double, double)> K, int nodes, int var) {
		this->a = a;
		this->b = b;
		this->f = f;
		this->K = K;
		this->nodes = nodes;
		this->var = var;
		switch (nodes)
		{
		case 5: // формула Боде
			G = { 28.0 / 45.0, 64.0 / 45.0, 24.0 / 45.0, 64.0 / 45.0, 14.0 / 45.0 };
			h = (b - a) / (nodes - 1);
			p = 7;
			break;
		case 4: // формула трех восьмых
			G = { 6.0 / 8.0, 9.0 / 8.0, 9.0 / 8.0, 3.0 / 8.0 };
			h = (b - a) / (nodes - 1);
			p = 5;
			break;
		case 3: // формула Симпсона
			G = { 2.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0 };
			h = (b - a) / (nodes - 1);
			p = 5;
			break;
		case 2: // формула трапеций
			G = { 2.0 / 2.0, 1.0 / 2.0 };
			h = (b - a) / (nodes - 1);
			p = 3;
			break;
		case 1: // формула средних прямоугольников
			G = { 2.0 / 2.0 };
			h = (b - a);
			p = 3;
			break;
		}
	}
	void solve();
	double computeIntegral();
	void PlotAdaptiveMesh();
	void Plot(int n);

};

