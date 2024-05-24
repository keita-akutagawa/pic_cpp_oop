#include <vector>
#include <cmath>
#include <iostream>

using namespace std;


void time_evolution_B(vector<vector<double>>& B, 
                      const vector<vector<double>> E, 
                      int n_x, double dx, double dt) 
{
    for (int i = 0; i < n_x; i++) {
            B[0][i] += 0.0;
            B[1][i] += (E[2][(i+1)%n_x] - E[2][i]) / dx * dt;
            B[2][i] += -(E[1][(i+1)%n_x] - E[1][i]) / dx * dt;
    }
}


void time_evolution_E(vector<vector<double>>& E, 
                      const vector<vector<double>> B, 
                      const vector<vector<double>> current,
                      int n_x, double dx, double dt, 
                      double c, double epsilon0) 
{
    for (int i = 0; i < n_x; i++) {
            E[0][i] += (-current[0][i] / epsilon0) * dt;
            E[1][i] += (-current[1][i] / epsilon0 
                        - pow(c, 2) * (B[2][i] - B[2][(n_x+i-1)%n_x]) / dx) * dt;
            E[2][i] += (-current[2][i] / epsilon0 
                        + pow(c, 2) * ((B[1][i] - B[1][(n_x+i-1)%n_x]) / dx)) * dt;
    }
}
