#include <vector>
#include <cmath>
#include <omp.h>
#include <iostream>

using namespace std;


void get_rho(vector<double>& rho, 
             vector<int>& r_index, 
             vector<double>& cr, 
             const vector<double> r, 
             int n_start, int n_last, double q, 
             int n_x, double dx)
{

    double cx1, cx2;
    int x_index_1, x_index_2;

    #pragma omp parallel for
    for (int i = 3 * n_start; i < 3 * n_last; i+=3) {
        r_index[i] = floor(r[i] / dx);
        cr[i] = r[i] / dx - r_index[i];
    }

    for (int i = 3 * n_start; i < 3 * n_last; i+=3) {

        cx1 = cr[i];
        cx2 = 1.0 - cr[i];

        x_index_1 = r_index[i];
        x_index_2 = x_index_1 + 1;

        if (x_index_2 != n_x) {

            rho[x_index_1] += q * cx2;
            rho[x_index_2] += q * cx1;

        } else {

            rho[x_index_1] += q * cx2;
            rho[0] += q * cx1;

        }
    }
    
}





