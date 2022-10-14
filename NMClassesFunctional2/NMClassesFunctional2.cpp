#include "FredholmEquationII.h"
#include "functions.h"
#include <excpt.h>

double f(double x) {
    return x - 0.5;
}
double K(double x, double t) {
    return 1;
}

int main(){
    std::function<double(double)> ft = f;
    std::function<double(double,double)> Kt = K;

    std::vector<std::function<double(double)>> F;
    std::vector<std::function<double(double, double)>> K;

    F.push_back(f11); K.push_back(K11);
    F.push_back(f12); K.push_back(K12);
    F.push_back(f13); K.push_back(K13);
    F.push_back(f14); K.push_back(K14);
    F.push_back(f15); K.push_back(K15);
    F.push_back(f16); K.push_back(K16);
    F.push_back(f17); K.push_back(K17);
    F.push_back(f18); K.push_back(K18);
    F.push_back(f19); K.push_back(K19);
    F.push_back(f20); K.push_back(K20);
    
    std::vector<FredholmEquationII> task;

    task.push_back(FredholmEquationII(0, 1, F[0], K[0], 4, 0));
    task.push_back(FredholmEquationII(1e-12, 0.5, F[1], K[1], 5, 1)); // wrong
    task.push_back(FredholmEquationII(0, 2 * PI , F[2], K[2], 5, 2));
    task.push_back(FredholmEquationII(0, 1, F[3], K[3], 5, 3)); 
    task.push_back(FredholmEquationII(0, 2 * PI, F[4], K[4], 5, 4));
    task.push_back(FredholmEquationII(0, 1, F[5], K[5], 4, 5)); // dont work
    task.push_back(FredholmEquationII(0, 1, F[6], K[6], 5, 6));
    task.push_back(FredholmEquationII(0, 1, F[7], K[7], 5, 7));
    task.push_back(FredholmEquationII(-1, 1, F[8], K[8], 5, 8)); // slow - small precision
    task.push_back(FredholmEquationII(0, PI, F[9], K[9], 5, 9)); // need small precision
    


    std::vector<FredholmEquationII> taskbynodes;

    taskbynodes.push_back(FredholmEquationII(0, 1, F[0], K[0], 5, 0));
    taskbynodes.push_back(FredholmEquationII(0, 1, F[0], K[0], 4, 0));
    taskbynodes.push_back(FredholmEquationII(0, 1, F[0], K[0], 3, 0));
    taskbynodes.push_back(FredholmEquationII(0, 1, F[0], K[0], 2, 0));
    taskbynodes.push_back(FredholmEquationII(0, 1, F[0], K[0], 1, 0));
    
    
    FredholmEquationII d(0, 1, ft, Kt, 5, 10);
    d.solve();
    d.computeIntegral();
    d.PlotAdaptiveMesh();
    //d.Plot(1024);
    std::system("cmd");
    return 0;
}
