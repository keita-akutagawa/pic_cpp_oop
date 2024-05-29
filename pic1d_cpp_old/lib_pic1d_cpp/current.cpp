#include <vector>
#include <cmath>
#include <omp.h>
#include <iostream>

using namespace std;


void get_current_density(vector<vector<double>>& current_tmp, 
                         vector<int>& r_index, 
                         vector<double>& cr, 
                         vector<double>& gamma,
                         const vector<double> r, 
                         const vector<double> v, 
                         int n_start, int n_last, 
                         double q, int n_x, double dx, double c)
{
    #pragma omp parallel for
    for (int i = 3 * n_start; i < 3 * n_last; i+=3) {
        r_index[i] = floor(r[i] / dx);
        cr[i] = r[i] / dx - r_index[i];
    }

    double tmp, tmp0, tmp1, tmp2;
    double cx1, cx2;
    int x_index_1, x_index_2;
    
    for (int i = 3 * n_start; i < 3 * n_last; i+=3) {

        tmp = q / gamma[i/3];
        tmp0 = tmp * v[i];
        tmp1 = tmp * v[i+1];
        tmp2 = tmp * v[i+2];

        cx1 = cr[i];
        cx2 = 1.0 - cr[i];

        x_index_1 = r_index[i];
        x_index_2 = x_index_1 + 1;

        if (x_index_2 != n_x) {

            current_tmp[0][x_index_1] += tmp0 * cx2;
            current_tmp[0][x_index_2] += tmp0 * cx1;

            current_tmp[1][x_index_1] += tmp1 * cx2;
            current_tmp[1][x_index_2] += tmp1 * cx1;

            current_tmp[2][x_index_1] += tmp2 * cx2;
            current_tmp[2][x_index_2] += tmp2 * cx1;

        } else {

            current_tmp[0][x_index_1] += tmp0 * cx2;
            current_tmp[0][0] += tmp0 * cx1;

            current_tmp[1][x_index_1] += tmp1 * cx2;
            current_tmp[1][0] += tmp1 * cx1;

            current_tmp[2][x_index_1] += tmp2 * cx2;
            current_tmp[2][0] += tmp2 * cx1;
        }
    }
}



