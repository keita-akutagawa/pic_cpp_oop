#include <vector>
#include <cmath>
#include <omp.h>
#include <algorithm>
#include <numeric>
#include <iostream>

using namespace std;


void sort_particles(vector<int>& index_species,
                    vector<double>& r, vector<double>& v,   
                    vector<double>& tmp_copy_double, 
                    int n_start, int n_last, 
                    int n_x, int n_y)
{
    #pragma omp parallel for
    for (int i = 3 * n_start; i < 3 * n_last; i+=3) {
        tmp_copy_double[i/3] = r[i];
    }

    //メモリ節約のため。x_indexにソートインデックスを保管。
    iota(index_species.begin(), index_species.end(), n_start);
    sort(index_species.begin(), index_species.end(), 
         [&](int a, int b){return tmp_copy_double[a] < tmp_copy_double[b];});

    //rのソート
    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {
        tmp_copy_double[i/3] = r[i];
    }
    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {
        r[i] = tmp_copy_double[index_species[i/3 - n_start]];
    }
    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {
        tmp_copy_double[i/3] = r[i+1];
    }
    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {
        r[i+1] = tmp_copy_double[index_species[i/3-n_start]];
    }
    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {
        tmp_copy_double[i/3] = r[i+2];
    }
    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {
        r[i+2] = tmp_copy_double[index_species[i/3-n_start]];
    }

    //vのソート
    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {
        tmp_copy_double[i/3] = v[i];
    }
    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {
        v[i] = tmp_copy_double[index_species[i/3 - n_start]];
    }
    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {
        tmp_copy_double[i/3] = v[i+1];
    }
    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {
        v[i+1] = tmp_copy_double[index_species[i/3-n_start]];
    }
    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {
        tmp_copy_double[i/3] = v[i+2];
    }
    #pragma omp parallel for
    for(int i = 3 * n_start; i < 3 * n_last; i+=3) {
        v[i+2] = tmp_copy_double[index_species[i/3-n_start]];
    }
}


