#include <vector>
#include <cmath>
#include <omp.h>

using namespace std;


void time_evolution_v(vector<double>& v, 
                      vector<double>& gamma, 
                      vector<double>& T, 
                      vector<double>& S, 
                      vector<double>& v_minus, 
                      vector<double>& v_0, 
                      vector<double>& v_plus, 
                      const vector<double> B_particle,
                      const vector<double> E_particle, 
                      int n_start, int n_last, 
                      double m, double q, double dt, double c)
{
    double tmp1 = q / m * (dt / 2.0);
    double tmp2 = 0.0;
    double tmp3 = 0.0;
    double tmp4 = 1.0 / (c * c);
    
    #pragma omp parallel for private(tmp2, tmp3)
    for (int i = 3 * n_start; i < 3 * n_last; i+=3) {
        gamma[i/3] = sqrt(1.0 + (v[i] * v[i] + v[i+1] * v[i+1] + v[i+2] * v[i+2]) * tmp4);

        tmp2 = tmp1 / gamma[i/3];
        T[i] = tmp2 * B_particle[i];
        T[i+1] = tmp2 * B_particle[i+1];
        T[i+2] = tmp2 * B_particle[i+2];
    
        tmp3 = 2.0 / (1.0 + T[i] * T[i] + T[i+1] * T[i+1] + T[i+2] * T[i+2]);
        S[i] = tmp3 * T[i];
        S[i+1] = tmp3 * T[i+1];
        S[i+2] = tmp3 * T[i+2];
    
        v_minus[i] = v[i] + tmp1 * E_particle[i];
        v_minus[i+1] = v[i+1] + tmp1 * E_particle[i+1];
        v_minus[i+2] = v[i+2] + tmp1 * E_particle[i+2];
    
        v_0[i] = v_minus[i] + (v_minus[i+1] * T[i+2] - v_minus[i+2] * T[i+1]);
        v_0[i+1] = v_minus[i+1] + (v_minus[i+2] * T[i] - v_minus[i] * T[i+2]);
        v_0[i+2] = v_minus[i+2] + (v_minus[i] * T[i+1] - v_minus[i+1] * T[i]);
    
        v_plus[i] = v_minus[i] + (v_0[i+1] * S[i+2] - v_0[i+2] * S[i+1]);
        v_plus[i+1] = v_minus[i+1] + (v_0[i+2] * S[i] - v_0[i] * S[i+2]);
        v_plus[i+2] = v_minus[i+2] + (v_0[i] * S[i+1] - v_0[i+1] * S[i]);
    
        v[i] = v_plus[i] + tmp1 * E_particle[i];
        v[i+1] = v_plus[i+1] + tmp1 * E_particle[i+1];
        v[i+2] = v_plus[i+2] + tmp1 * E_particle[i+2];
    }
}


void time_evolution_x(vector<double>& r, 
                      vector<double>& gamma,
                      const vector<double> v,
                      int n_particle, double dt, double c)
{
    #pragma omp parallel for
    for (int i = 0; i < 3 * n_particle; i+=3){
        double tmp = dt / gamma[i/3];
        r[i] += tmp * v[i];
        r[i+1] += tmp * v[i+1];
        r[i+2] += tmp * v[i+2];
    }
}

