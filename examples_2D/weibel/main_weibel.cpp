#include <iostream>
#include <iomanip>
#include <random>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <omp.h>
//#include <filesystem>
#include <numeric>
#include <chrono>
#include "module.hpp"

using namespace std;


int main()
{
    string dirname = "weibel", filename = "weibel";

    const double c = 0.5, epsilon0 = 1.0, mu_0 = 1.0 / (epsilon0 * pow(c, 2));
    const int n_x = 256, n_y = 256, n_e = 20, n_i = 20, step = 3000;
    const int n_electron = n_e * n_x * n_y, n_ion = n_i * n_x * n_y;
    const int n_particle = n_electron + n_ion;
    const int iteration = 10000;
    double dx = 1.0, dy = 1.0, x_max = n_x * dx, y_max = n_y * dy, dt = 1.0;

    //if(!filesystem::exists(dirname)) {
    //    filesystem::create_directory(dirname);
    //}

    cout << "Total number of particles is " << n_particle << ".\n";

    double B0;
    double t_r, T_e, T_i;
    double m_unit, r_m, m_electron, m_ion;
    double q_unit, r_q, q_ion, q_electron;
    double omega_pe, omega_pi, omega_ce, omega_ci;
    double V_A, C_S, debye_length;
    vector<double> v_ion(3, 0), v_electron(3, 0);
    double v_thi, v_the;

    
    m_unit = 1.0;
    r_m = 1.0;
    t_r = 1.0;
    m_electron = 1 * m_unit;
    m_ion = m_electron / r_m;
    r_q = 1.0;
    T_e = 1.0/2.0 * m_electron * pow(0.1*c, 2);
    T_i = 1.0/2.0 * m_ion * pow(0.1*c, 2);
    B0 = 1.0;
    q_unit = sqrt(epsilon0 * T_e / static_cast<double>(n_e));
    q_electron = -1 * q_unit;
    q_ion = r_q * q_unit;
    omega_pe = sqrt(static_cast<double>(n_e) * pow(q_electron, 2) / m_electron / epsilon0);
    omega_pi = sqrt(static_cast<double>(n_i) * pow(q_ion, 2) / m_ion / epsilon0);
    omega_ce = q_electron * B0 / m_electron;
    omega_ci = q_ion * B0 / m_ion;
    V_A = B0 / sqrt(mu_0 * (static_cast<double>(n_i) * m_ion));
    C_S = sqrt(r_m * T_e);
    debye_length = sqrt(epsilon0 * T_e / static_cast<double>(n_e) / pow(q_electron, 2));
    dx = debye_length;
    dy = debye_length;
    x_max = n_x * dx;
    y_max = n_y * dy;
    v_thi = sqrt(T_i / m_ion);
    v_the = sqrt(T_e / m_electron);
    v_ion[0] = 0.0;
    v_ion[1] = 0.0;
    v_ion[2] = 0.0;
    v_electron[0] = 0.0;
    v_electron[1] = 0.0;
    v_electron[2] = 0.0;
    
///////////////////////////////////////////////////////////////////////

    vector<vector<vector<double>>> B(3, vector<vector<double>>(n_x, vector<double>(n_y, 0.0)));
    vector<vector<vector<double>>> B_tmp(3, vector<vector<double>>(n_x, vector<double>(n_y, 0.0)));
    vector<vector<vector<double>>> E(3, vector<vector<double>>(n_x, vector<double>(n_y, 0.0)));
    vector<vector<vector<double>>> E_tmp(3, vector<vector<double>>(n_x, vector<double>(n_y, 0.0)));
    vector<vector<vector<double>>> current(3, vector<vector<double>>(n_x, vector<double>(n_y, 0.0)));
    vector<vector<vector<double>>> current_tmp(3, vector<vector<double>>(n_x, vector<double>(n_y, 0.0)));
    vector<double> r(3 * n_particle, 0.0);
    vector<double> v(3 * n_particle, 0.0);
    vector<double> gamma(n_particle, 0.0);
    vector<double> B_particle(3 * n_particle, 0.0);
    vector<double> E_particle(3 * n_particle, 0.0);
    vector<vector<double>> rho(n_x, vector<double>(n_y, 0.0));
    vector<vector<double>> phi(n_x, vector<double>(n_y, 0.0));
    vector<int> index_ion(n_ion, 0.0);
    vector<int> index_electron(n_electron, 0.0);
    vector<int> r_index(3 * n_particle, 0);
    vector<double> cr(3 * n_particle, 0.0);
    vector<double> v_minus(3 * n_particle, 0.0);
    vector<double> v_0(3 * n_particle, 0.0);
    vector<double> v_plus(3 * n_particle, 0.0);
    vector<double> S(3 * n_particle, 0.0);
    vector<double> T(3 * n_particle, 0.0);
    vector<double> tmp_copy_double(n_particle, 0.0);

    ///////////////////////////////////////////////////////////////////////

    set_initial_position_x(r, 0, n_ion, x_max, 1);
    set_initial_position_y(r, 0, n_ion, y_max, 2);
    set_initial_position_x(r, n_ion, n_ion + n_electron, x_max, 4);
    set_initial_position_y(r, n_ion, n_ion + n_electron, y_max, 5);

    set_initial_velocity_x(v, 0, n_ion, v_ion[0], v_thi, 7); 
    set_initial_velocity_y(v, 0, n_ion, v_ion[1], v_thi, 8); 
    set_initial_velocity_z(v, 0, n_ion, v_ion[2], v_thi*5.0, 9); 
    set_initial_velocity_x(v, n_ion, n_ion + n_electron, v_electron[0], v_the, 10); 
    set_initial_velocity_y(v, n_ion, n_ion + n_electron, v_electron[1], v_the, 11); 
    set_initial_velocity_z(v, n_ion, n_ion + n_electron, v_electron[2], v_the*5.0, 12);  

    sort_particles(index_ion, r, v, 
                   tmp_copy_double, 0, n_ion, n_x, n_y);
    sort_particles(index_electron, r, v, 
                   tmp_copy_double, n_ion, n_ion+n_electron, n_x, n_y);

    get_rho(rho, r_index, cr, r, 
            0, n_ion, q_ion, n_x, n_y, dx, dy);
    get_rho(rho, r_index, cr, r, 
            n_ion, n_ion+n_electron, q_electron, n_x, n_y, dx, dy);
    poisson_solver_jacobi(phi, rho, iteration, n_x, n_y, dx, dy, epsilon0);
    get_E(E, phi, n_x, n_y, dx, dy);


    chrono::system_clock::time_point  start, end, start_total, end_total;
    double elapsed;
    string path;
    path = "time.txt";
    ofstream file_time(path);

    for (int k = 0; k < step+1; k++) {
        start_total = chrono::system_clock::now();

        if (k % 10 == 0) {
            cout << int(k*dt) << "step done..." << "\n";
        }

        //OUTPUT
        if (k % 100 == 0) {
            io_xvKE(r, v, 
                    dirname, filename, k, 
                    n_ion, n_electron, n_x, n_y, m_ion, m_electron);
            io_EBcurrent(E, B, current, 
                  dirname, filename, k, 
                  n_x, n_y);
        }     

        //STEP2
        time_evolution_B(B, E, n_x, n_y, dx, dy, dt/2.0);

        //STEP3
        //電場・磁場を整数格子点上に再定義
        for (int i = 0; i < n_x; i++) {
            for (int j = 0; j < n_y; j++) {
                B_tmp[0][i][j] = (B[0][i][j] + B[0][i][(n_y+j-1)%n_y]) / 2.0;
                B_tmp[1][i][j] = (B[1][i][j] + B[1][(n_x+i-1)%n_x][j]) / 2.0;
                B_tmp[2][i][j] = (B[2][i][j] + B[2][(n_x+i-1)%n_x][j] + B[2][i][(n_y+j-1)%n_y] + B[2][(n_x+i-1)%n_x][(n_y+j-1)%n_y]) / 4.0;
                E_tmp[0][i][j] = (E[0][i][j] + E[0][(n_x+i-1)%n_x][j]) / 2.0;
                E_tmp[1][i][j] = (E[1][i][j] + E[1][i][(n_y+j-1)%n_y]) / 2.0;
                E_tmp[2][i][j] = E[2][i][j];
            }
        }
        //start = chrono::system_clock::now();
        get_particle_field(B_particle, E_particle, 
                           r_index, cr,  
                           B_tmp, E_tmp, r, 0, n_particle, n_x, n_y, dx, dy);
        //end = chrono::system_clock::now();
        //elapsed = chrono::duration_cast<chrono::milliseconds>(end-start).count();
        //cout << "field interpolation time is " << elapsed/1000.0 << "\n";
        //start = chrono::system_clock::now();
        time_evolution_v(v, gamma, T, S, v_minus, v_0, v_plus, 
                         B_particle, E_particle, 
                         0, n_ion, 
                         m_ion, q_ion, dt, c);
        time_evolution_v(v, gamma, T, S, v_minus, v_0, v_plus, 
                         B_particle, E_particle, 
                         n_ion, n_ion + n_electron, 
                         m_electron, q_electron, dt, c);
        //end = chrono::system_clock::now();
        //elapsed = chrono::duration_cast<chrono::milliseconds>(end-start).count();
        //cout << "velocity time is " << elapsed/1000.0 << "\n";

        //STEP4
        time_evolution_x(r, gamma, v, n_particle, dt/2.0, c);
        periodic_boudary_condition_x(r, n_particle, x_max);
        periodic_boudary_condition_y(r, n_particle, y_max);

        //STEP5
        //0で初期化
        for (int comp = 0; comp < 3; comp++) {
            for (int i = 0; i < n_x; i++) {
                for (int j = 0; j < n_y; j++) {
                    current_tmp[comp][i][j] = 0.0;
                }
            }
        }
        //start = chrono::system_clock::now();
        get_current_density(current_tmp, r_index, 
                            cr, gamma, r, v, 
                            0, n_ion, q_ion, n_x, n_y, dx, dy, c);
        get_current_density(current_tmp, r_index, 
                            cr, gamma, r, v, 
                            n_ion, n_ion + n_electron, q_electron, n_x, n_y, dx, dy, c);
        //end = chrono::system_clock::now();
        //elapsed = chrono::duration_cast<chrono::milliseconds>(end-start).count();
        //cout << "current density time is " << elapsed/1000.0 << "\n";
        //半整数格子点上に再定義
        for (int i = 0; i < n_x; i++) {
            for (int j = 0; j < n_y; j++) {
                current[0][i][j] = (current_tmp[0][i][j] + current_tmp[0][(i+1)%n_x][j]) / 2.0;
                current[1][i][j] = (current_tmp[1][i][j] + current_tmp[1][i][(j+1)%n_y]) / 2.0;
                current[2][i][j] = current_tmp[2][i][j];
            }
        }

        //STEP6
        time_evolution_B(B, E, n_x, n_y, dx, dy, dt/2.0);
        
        //STEP7
        time_evolution_E(E, B, current, n_x, n_y, dx, dy, dt, c, epsilon0);
        
        //STEP8
        time_evolution_x(r, gamma, v, n_particle, dt/2.0, c);
        periodic_boudary_condition_x(r, n_particle, x_max);
        periodic_boudary_condition_y(r, n_particle, y_max);

        if (k % 50 == 0) {
            sort_particles(index_ion, r, v, 
                   tmp_copy_double, 0, n_ion, n_x, n_y);
            sort_particles(index_electron, r, v, 
                        tmp_copy_double, n_ion, n_ion+n_electron, n_x, n_y);
        } 

        end_total = chrono::system_clock::now();
        elapsed = chrono::duration_cast<chrono::milliseconds>(end_total-start_total).count();
        cout << "total time is " << elapsed / 1000.0 << "\n";
        cout << "------------------------------\n";

        file_time << elapsed / 1000.0 << "\n";

    }
}