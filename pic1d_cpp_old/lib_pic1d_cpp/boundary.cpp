#include <vector>
#include <cmath>
#include <omp.h>

using namespace std;


void periodic_boudary_condition_x(vector<double>& r,  
                                  int n_particle, double x_max)
{
    #pragma omp parallel for
    for (int i = 0; i < 3 * n_particle; i+=3){
        if (r[i] > x_max) {
            r[i] = 1e-10;
        } else if (r[i] < 0.0) {
            r[i] = x_max - 1e-10;
        }
    }
}



void refrective_boudary_condition_x(vector<double>& v, vector<double>& r,  
                                    int n_particle, double x_min, double x_max)
{
    #pragma omp parallel for
    for (int i = 0; i < 3 * n_particle; i+=3){
        if (r[i] > x_max) {
            v[i] = -v[i];
            r[i] = x_max - 1e-10;
        } else if (r[i] < x_min) {
            v[i] = -v[i];
            r[i] = x_min + 1e-10;
        }
    }
}


