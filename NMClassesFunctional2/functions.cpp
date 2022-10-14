#include"functions.h"

double f11(double x) { return exp(-x); }
double K11(double x, double t) { return x * exp(t) / 2; }

double f12(double x) { return 1 + (cos(x / 2) - 1) / x; }
double K12(double x, double t) { return sin(x*t); }

double f13(double x) { return (3 * cos(2 * x) + 5) / (16*3.14); }
double K13(double x, double t) { return -1.0 / (4*3.14*(pow(sin((x+t)/2),2) + 0.25 * pow(cos((x + t) / 2), 2))); }

double f14(double x) { return x / 6; }
double K14(double x, double t) { return 2 * lambda * x - t; }

double f15(double x) { return cos(2*x); }
double K15(double x, double t) { return sin(x) * cos(t); }

double f16(double x) { return x; }
double K16(double x, double t) { return lambda * (4 * x * t - x * x); }

double f17(double x) { return 1; }
double K17(double x, double t) { return x * t * t; }

double f18(double x) { return 5*x/6; }
double K18(double x, double t) { return (x * t)/2; }

double f19(double x) { return 1 - x * (exp(x)-exp(-x)); }
double K19(double x, double t) { return x * x * exp(x * t); }

double f20(double x) { return 1; }
double K20(double x, double t) { return lambda * cos(x)*cos(x); }
