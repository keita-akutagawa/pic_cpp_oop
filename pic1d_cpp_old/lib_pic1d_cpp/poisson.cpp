#include <vector>
#include <cmath>
#include <omp.h>
#include <iostream>

using namespace std;


void poisson_solver_jacobi(vector<double>& phi, 
                           const vector<double> rho, int iteration,
                           int n_x, double dx, double epsilon0)
{
    double tmp = 1.0 / (2.0 * (1.0/pow(dx, 2)));
    double tmp_x = 1.0 / pow(dx, 2);
    for (int iter = 0; iter < iteration; iter++) {
        for (int i = 1; i < n_x-1; i++) {
                phi[i] = tmp * (rho[i]/epsilon0 
                          + (phi[i-1] + phi[i+1]) * tmp_x);
        }
    }
}


void get_E(vector<vector<double>>& E, vector<double> phi, 
           int n_x, double dx)
{
    for (int i = 0; i < n_x; i++) {
            E[0][i] = -(phi[(i+1)%n_x] - phi[i]) / dx;
    }
}