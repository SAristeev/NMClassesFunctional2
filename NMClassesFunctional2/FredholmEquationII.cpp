#include "FredholmEquationII.h"
constexpr auto ABSEPS = 1e-12;
constexpr auto RELEPS = 1e-12;

auto max = [](auto x, auto y) {return x > y ? x : y; };

void LUdecomposition(std::vector<std::vector<double>> x,
	std::vector<std::vector<double>>& lower,
	std::vector<std::vector<double>>& upper) {
	size_t n = x.size();
	for (int i = 0; i < n; i++) {
		std::vector<double> tmpl, tmpu;
		for (int j = 0; j < n; j++) {
			tmpl.push_back(i == j ? 1 : 0);
			tmpu.push_back(0);
		}
		lower.push_back(tmpl);
		upper.push_back(tmpu);
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i <= j) {
				double sum = 0;
				for (int k = 0; k < i; k++) {
					sum += lower[i][k] * upper[k][j];
				}
				upper[i][j] = x[i][j] - sum;
			}
			if (i > j) {
				double sum = 0;
				for (int k = 0; k < j; k++) {
					sum += lower[i][k] * upper[k][j];
				}
				lower[i][j] = (x[i][j] - sum) / upper[j][j];
			}
		}
	}
}

std::vector<double> SolveLinearSystemLU(std::vector<std::vector<double>> A, std::vector<double> b) {
	size_t n = A.size();
	std::vector<std::vector<double>> L, U;
	LUdecomposition(A, L, U);
	std::vector<double> x(n), y(n);
	for (int i = 0; i < n; i++) {
		double sum = 0;
		for (int j = 0; j < i; j++) {
			sum += L[i][j] * y[j];
		}
		y[i] = b[i] - sum;
	}
	for (int i = n - 1; i >= 0; --i) {
		double sum = 0;
		for (int j = i + 1; j < n; j++) {
			sum += U[i] [j] * x[j];
		}
		x[i] = (y[i] - sum) / U[i][i];
	}
	return x;
}

double AdaptiveIntegrate(std::function<double(double)> u, double a, double b, std::vector<double> G, int nodes, int p, std::set<double> &mesh) {
	double alpha = a, beta = b;
	double delta = 0, hlocal = (beta - alpha) / (nodes - 1);
	double I = 0, Ih = 0, Ih2 = 0, E = 0, deep = 0, deepmax = 5;
	mesh.insert(alpha);
	mesh.insert(beta);
	auto LocalI = [u, G, nodes](double a, double b) {
		double h = nodes != 1 ? ((b - a) / (nodes - 1)) : ((b - a) / 2);
		double tmp = G[nodes - 1] * u(a) * h;
		for (int i = 1; i < nodes; i++) {
			tmp += G[i] * u(a + i * h) * h;
		}
		return tmp;
	};
	while (beta - alpha > ABSEPS)
	{
		Ih = LocalI(alpha, beta);
		Ih2 = LocalI(alpha, beta/2) + LocalI(beta / 2, beta);
		delta = (Ih2 - Ih) / (pow(2, p) - 1);
		if (abs(delta) < max(ABSEPS, RELEPS * abs(I)) * (beta - alpha) / (nodes - 1) * (b - a)) {
			deep = 0;
			I += Ih2 + delta;
			E += delta;
			alpha = beta;
			beta = b;
		}
		else {
			beta = (alpha + beta) / 2;
			mesh.insert(beta);
			deep++;
			if (deep > deepmax) {
				deep = 0;
				I += Ih2 + delta;
				E += delta;
				alpha = beta;
				beta = b;
			}
		}
	}
	return I;
}

double AdaptiveIntegrate(std::function<double(double)> u, double a, double b, std::vector<double> G, int nodes, int p) {
	double alpha = a, beta = b;
	double delta = 0, hlocal = (beta - alpha) / (nodes - 1);
	double I = 0, Ih = 0, Ih2 = 0, E = 0, deep = 0, deepmax = 5;
	auto LocalI = [u, G, nodes](double a, double b) {
		double h = nodes != 1 ? ((b - a) / (nodes - 1)) : ((b - a) / 2);
		double tmp = G[nodes - 1] * u(a) * h;
		for (int i = 1; i < nodes; i++) {
			tmp += G[i] * u(a + i * h) * h;
		}
		return tmp;
	};
	while (beta - alpha > ABSEPS)
	{
		Ih = LocalI(alpha, beta);
		Ih2 = LocalI(alpha, beta / 2) + LocalI(beta / 2, beta);
		delta = (Ih2 - Ih) / (pow(2, p) - 1);
		if (abs(delta) < max(ABSEPS, RELEPS * abs(I)) * (beta - alpha) / (nodes - 1) * (b - a)) {
			deep = 0;
			I += Ih2 + delta;
			E += delta;
			alpha = beta;
			beta = b;
		}
		else {
			beta = (alpha + beta) / 2;
			deep++;
			if (deep > deepmax) {
				deep = 0;
				I += Ih2 + delta;
				E += delta;
				alpha = beta;
				beta = b;
			}
		}
	}
	return I;
}

void FredholmEquationII::solve(){
	int M = nodes;
	std::function<double(double)> uprev = std::function<double(double)>([](double x) { return 1; });
	double hM;
	while (abs(AdaptiveIntegrate(std::function<double(double)>([u = u, uprev](double x) {return u(x) - uprev(x); }), a, b, G, nodes, p)) > 1e-2) {
		M *= 2;
		hM = (b - a) / (M - 1);
		uprev = u;

		auto index = [nodes = nodes, M](int i) {
			if (nodes == 1) {
				return 0;
			}
			else {
				if (i == 0 || i == M - 1) {
					return nodes - 1;
				}
				else
					return i % (nodes - 1);
			}
		};
		
		std::vector<std::vector<double>> A;

		for (int i = 0; i < M; i++) {
			std::vector<double> column;
			for (int j = 0; j < M; j++) {
				double tmp = i == j ? 1 : 0;
				tmp -= G[index(j)] * K(a + i * hM, a + j * hM) * hM;
				column.push_back(tmp);
			}
			A.push_back(column);
		}
		std::ofstream File;
		std::string FileName = "A" + std::to_string(M) + ".txt";
		File.open(FileName);
		for (int i = 0; i < M; i++) {			
			for (int j = 0; j < M; j++) {
				File << std::setw(10) << A[i][j] << " ";
			}
			File << std::endl;
		}
		File.close();
		std::cout << std::endl;
		std::vector<double> F;
		for (int i = 0; i < M; i++) {
			F.push_back(f(a + i * hM));
		}

		U = SolveLinearSystemLU(A, F);

		u = std::function<double(double)>([a = a, b = b, f = f, K = K, U=U, G=G, nodes = nodes, index = index, M = M](double x) {
			if (U.size() == 0) {
				return 0.0;
			}
			double tmp = 0;
			double h = (b - a) / (M - 1);
			for (int i = 0; i < M; i++) {
				tmp += G[index(i)] * K(x, a + i * h) * U[i] * h;
			}
			return tmp + f(x);
			}
		);
	}
}

double FredholmEquationII::computeIntegral() {
	I = AdaptiveIntegrate(u, a, b, G, nodes, p, mesh);
	return I;
}

void FredholmEquationII::PlotAdaptiveMesh() {
	std::ofstream FileX, FileU, FileV;

	FileV.open("V.txt");
	FileX.open("X.txt");
	FileU.open("U.txt");

	for (auto x : mesh) {
		FileX << x << ", ";
		FileU << u(x) << ", ";
	}

	FileX << std::endl;
	FileU << std::endl;

	FileV << var + 1 << " " << lambda << " " << I << " " << mesh.size();

	FileX.close();
	FileU.close();
	FileV.close();
	std::system("python plot.py");
};
void FredholmEquationII::Plot(int n) {
	std::ofstream FileX, FileU, FileV;
	FileX.open("X.txt");
	FileU.open("U.txt");
	FileV.open("V.txt");

	for (int i = 0; i < n; i++) {
		FileX << a + i * (b - a) / (n - 1) << ", ";
		FileU << u(a + i * (b - a) / (n - 1)) << ", ";
	}
	FileX << std::endl;
	FileU << std::endl;
	FileV << var + 1 << " " << lambda << " " << I << " " << mesh.size();


	FileX.close();
	FileU.close();
	FileV.close();
	std::system("python plot.py");
};