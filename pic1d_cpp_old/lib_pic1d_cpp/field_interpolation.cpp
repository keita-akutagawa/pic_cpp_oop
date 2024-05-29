#include <vector>
#include <cmath>
#include <omp.h>
#include <iostream>
#include <algorithm>

using namespace std;


void get_particle_field(vector<double>& B_particle, 
                        vector<double>& E_particle, 
                        vector<int>& r_index, 
                        vector<double>& cr,  
                        const vector<vector<double>> B_tmp, 
                        const vector<vector<double>> E_tmp, 
                        const vector<double> r, 
                        int n_start, int n_last, 
                        int n_x, double dx)
{
    //加えていくので、先にリセットする。これを忘れると堆積していって発散する。
    fill(B_particle.begin(), B_particle.end(), 0.0);
    fill(E_particle.begin(), E_particle.end(), 0.0);


    double cx1, cx2;
    int x_index_1, x_index_2;

    #pragma omp parallel for private(cx1, cx2, x_index_1, x_index_2)
    for (int i = 3 * n_start; i < 3 * n_last; i+=3) {

        r_index[i] = floor(r[i] / dx);
        cr[i] = r[i] / dx - r_index[i];

        cx1 = cr[i];
        cx2 = 1.0 - cr[i];

        x_index_1 = r_index[i];
        x_index_2 = x_index_1 + 1;

        if (x_index_2 != n_x) {

            B_particle[i] += B_tmp[0][x_index_1] * cx2;
            B_particle[i] += B_tmp[0][x_index_2] * cx1;

            B_particle[i+1] += B_tmp[1][x_index_1] * cx2;
            B_particle[i+1] += B_tmp[1][x_index_2] * cx1;

            B_particle[i+2] += B_tmp[2][x_index_1] * cx2;
            B_particle[i+2] += B_tmp[2][x_index_2] * cx1;

            E_particle[i] += E_tmp[0][x_index_1] * cx2;
            E_particle[i] += E_tmp[0][x_index_2] * cx1;

            E_particle[i+1] += E_tmp[1][x_index_1] * cx2;
            E_particle[i+1] += E_tmp[1][x_index_2] * cx1;

            E_particle[i+2] += E_tmp[2][x_index_1] * cx2;
            E_particle[i+2] += E_tmp[2][x_index_2] * cx1;

        } else {

            B_particle[i] += B_tmp[0][x_index_1] * cx2;
            B_particle[i] += B_tmp[0][0] * cx1;

            B_particle[i+1] += B_tmp[1][x_index_1] * cx2;
            B_particle[i+1] += B_tmp[1][0] * cx1;

            B_particle[i+2] += B_tmp[2][x_index_1] * cx2;
            B_particle[i+2] += B_tmp[2][0] * cx1;

            E_particle[i] += E_tmp[0][x_index_1] * cx2;
            E_particle[i] += E_tmp[0][0] * cx1;
            
            E_particle[i+1] += E_tmp[1][x_index_1] * cx2;
            E_particle[i+1] += E_tmp[1][0] * cx1;

            E_particle[i+2] += E_tmp[2][x_index_1] * cx2;
            E_particle[i+2] += E_tmp[2][0] * cx1;
        }
    }
}

