#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;


void io_xvKE(const vector<double> r, 
             const vector<double> v, 
             string dirname, string filename, int k, 
             int n_ion, int n_electron, int n_x, double m_ion, double m_electron)
{
    //ios::sync_with_stdio(false);
    //cin.tie(nullptr);

    string path_x, path_v, path_KE;

    path_x = "./" + dirname + "/" + filename + "_x_" + to_string(k) + ".csv";
    path_v = "./" + dirname + "/" + filename + "_v_" + to_string(k) + ".csv";
    path_KE = "./" + dirname + "/" + filename + "_KE_" + to_string(k) + ".csv";
    ofstream file_x(path_x);
    ofstream file_v(path_v);
    ofstream file_KE(path_KE);
    double KE = 0.0;

    for (int i = 0; i < 3 * (n_ion + n_electron); i+=3) {
        file_x << r[i] << ',' << r[i+1] << ',' << r[i+2] << "\n";
        file_v << v[i] << ',' << v[i+1] << ',' << v[i+2] << "\n";

        if (i < 3 * n_ion) {
            KE += 1.0 / 2.0 * m_ion * (v[i] * v[i] + v[i+1] * v[i+1] + v[i+2] * v[i+2]);
        } else {
            KE += 1.0 / 2.0 * m_electron * (v[i] * v[i] + v[i+1] * v[i+1] + v[i+2] * v[i+2]);
        }
    }
    file_KE << KE;
}


void io_EBcurrent(const vector<vector<double>> E, 
                  const vector<vector<double>> B, 
                  const vector<vector<double>> current, 
                  string dirname, string filename, int k, 
                  int n_x)
{
    //ios::sync_with_stdio(false);
    //cin.tie(nullptr);

    string path_E, path_B, path_current;

    path_E = "./" + dirname + "/" + filename + "_E_" + to_string(k) + ".csv";
    path_B = "./" + dirname + "/" + filename + "_B_" + to_string(k) + ".csv";
    path_current = "./" + dirname + "/" + filename + "_current_" + to_string(k) + ".csv";
    ofstream file_E(path_E);
    ofstream file_B(path_B);
    ofstream file_current(path_current);
    
    for (int i = 0; i < n_x; i++) {
            file_E << E[0][i] << ',' << E[1][i] << ',' << E[2][i] << ',' 
            << i << ',' << 0 << ',' << 0 << "\n";
        
            file_B << B[0][i] << ',' << B[1][i] << ',' << B[2][i] << ','
            << i << ',' << 0 << ',' << 0 << "\n";

            file_current << current[0][i] << ',' << current[1][i] << ',' << current[2][i] << ','
            << i << ',' << 0 << ',' << 0 << "\n";
    }
}

