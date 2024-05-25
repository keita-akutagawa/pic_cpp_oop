#include <iostream>
#include <iomanip>
#include <random>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <omp.h>
#include <chrono>
#include "module.hpp"

using namespace std;


int main()
{
    string dirname = "two_stream_electron", filename = "two_stream_electron";

    const double c = 0.5, epsilon0 = 1.0, mu_0 = 1.0 / (epsilon0 * pow(c, 2));
    const int n_x = 512, n_e = 100, n_i = 100, step = 30000;
    const int n_electron = n_e * n_x / 2, n_beam = n_e * n_x / 2, n_ion = n_i * n_x;
    const int n_particle = n_electron + n_beam + n_ion;
    const int iteration = 10000;
    double dx = 1.0, x_max = n_x * dx, dt = 1.0;

    cout << "Total number of particles is " << n_particle << ".\n";

    double B0;
    double t_r, T_e, T_i;
    double m_unit, r_m, m_electron, m_ion;
    double q_unit, r_q, q_ion, q_electron;
    double omega_pe, omega_pi, omega_ce, omega_ci;
    double V_A, C_S, debye_length;
    vector<double> v_ion(3, 0), v_electron(3, 0), v_beam(3, 0);
    double v_thi, v_the;

    
    m_unit = 1.0;
    r_m = 0.01;
    t_r = 0.01;
    m_electron = 1 * m_unit;
    m_ion = m_electron / r_m;
    r_q = 1.0;
    T_e = 1.0/2.0 * m_electron * pow(0.01*c, 2);
    T_i = T_e / t_r;
    B0 = sqrt(static_cast<double>(n_e)) / 10.0;
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
    x_max = n_x * dx;
    v_thi = sqrt(T_i / m_ion);
    v_the = sqrt(T_e / m_electron);
    v_ion[0] = 0.0;
    v_ion[1] = 0.0;
    v_ion[2] = 0.0;
    v_electron[0] = -10.0 * v_thi;
    v_electron[1] = 0.0;
    v_electron[2] = 0.0;
    v_beam[0] = 10.0 * v_thi;
    v_beam[1] = 0.0;
    v_beam[2] = 0.0;
    
///////////////////////////////////////////////////////////////////////

    vector<vector<double>> B(3, vector<double>(n_x, 0.0));
    vector<vector<double>> B_tmp(3, vector<double>(n_x, 0.0));
    vector<vector<double>> E(3, vector<double>(n_x, 0.0));
    vector<vector<double>> E_tmp(3, vector<double>(n_x, 0.0));
    vector<vector<double>> current(3, vector<double>(n_x, 0.0));
    vector<vector<double>> current_tmp(3, vector<double>(n_x, 0.0));
    vector<double> r(3 * n_particle, 0.0);
    vector<double> v(3 * n_particle, 0.0);
    vector<double> gamma(n_particle, 0.0);
    vector<double> B_particle(3 * n_particle, 0.0);
    vector<double> E_particle(3 * n_particle, 0.0);
    vector<double> rho(n_x, 0.0);
    vector<double> phi(n_x, 0.0);
    vector<int> index_ion(n_ion, 0.0);
    vector<int> index_electron(n_electron+n_beam, 0.0);
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
    set_initial_position_x(r, n_ion, n_ion + n_electron + n_beam, x_max, 2);

    set_initial_velocity_x(v, 0, n_ion, v_ion[0], v_thi, 100); 
    set_initial_velocity_y(v, 0, n_ion, v_ion[1], v_thi, 101); 
    set_initial_velocity_z(v, 0, n_ion, v_ion[2], v_thi, 102); 
    set_initial_velocity_x(v, n_ion, n_ion + n_electron, v_electron[0], v_the, 103); 
    set_initial_velocity_y(v, n_ion, n_ion + n_electron, v_electron[1], v_the, 104); 
    set_initial_velocity_z(v, n_ion, n_ion + n_electron, v_electron[2], v_the, 105); 
    set_initial_velocity_x(v, n_ion + n_electron, n_ion + n_electron + n_beam, v_beam[0], v_the, 106);
    set_initial_velocity_y(v, n_ion + n_electron, n_ion + n_electron + n_beam, v_beam[1], v_the, 107); 
    set_initial_velocity_z(v, n_ion + n_electron, n_ion + n_electron + n_beam, v_beam[2], v_the, 108);  

    get_rho(rho, r_index, cr, r, 0, n_ion, q_ion, n_x, dx);
    get_rho(rho, r_index, cr, r, n_ion, n_ion + n_electron + n_beam, q_electron, n_x, dx);
    poisson_solver_jacobi(phi, rho, iteration, n_x, dx, epsilon0);
    get_E(E, phi, n_x, dx);


    for (int k = 0; k < step+1; k++) {
        if (k % 100 == 0) {
            cout << int(k*dt) << "step done..." << "\n";
        }
        //OUTPUT
        if (k % 100 == 0) {
            io_xvKE(r, v, 
               dirname, filename, k, 
               n_ion, n_electron + n_beam, n_x, m_ion, m_electron);
            io_EBcurrent(E, B, current, 
                         dirname, filename, k, 
                         n_x);
        }
        //STEP2
        time_evolution_B(B, E, n_x, dx, dt/2.0);
        //STEP3
        //電場・磁場を整数格子点上に再定義
        for (int i = 0; i < n_x; i++) {
            if (i == 0) {
                B_tmp[0][i] = B[0][i];
                B_tmp[1][i] = (B[1][i] + B[1][n_x-1]) / 2.0;
                B_tmp[2][i] = (B[2][i] + B[2][n_x-1]) / 2.0;
                E_tmp[0][i] = (E[0][i] + E[0][n_x-1]) / 2.0;
                E_tmp[1][i] = E[1][i];
                E_tmp[2][i] = E[2][i];
            } else {
            B_tmp[0][i] = B[0][i];
            B_tmp[1][i] = (B[1][i] + B[1][i-1]) / 2.0;
            B_tmp[2][i] = (B[2][i] + B[2][i-1]) / 2.0;
            E_tmp[0][i] = (E[0][i] + E[0][i-1]) / 2.0;
            E_tmp[1][i] = E[1][i];
            E_tmp[2][i] = E[2][i];
            }
        }
        get_particle_field(B_particle, E_particle, 
                           r_index, cr, 
                           B_tmp, E_tmp, r, 0, n_particle, n_x, dx);
        time_evolution_v(v, gamma, 
                         T, S, v_minus, v_0, v_plus, 
                         B_particle, E_particle, 
                         0, n_ion, 
                         m_ion, q_ion, dt, c);
        time_evolution_v(v, gamma, 
                         T, S, v_minus, v_0, v_plus, 
                         B_particle, E_particle, 
                         n_ion, n_ion + n_electron + n_beam, 
                         m_electron, q_electron, dt, c);
        //STEP4
        time_evolution_x(r, gamma, v, 
                         n_particle, dt/2.0, c);
        periodic_boudary_condition_x(r, n_particle, x_max);
        //STEP5
        //0で初期化
        for (int i = 0; i < 3; i++) {
            fill(current_tmp[i].begin(), current_tmp[i].end(), 0.0);
        }
        get_current_density(current_tmp, r_index, cr, gamma,
                            r, v, 
                            0, n_ion, q_ion, n_x, dx, c);
        get_current_density(current_tmp, r_index, cr, gamma, 
                            r, v,  
                            n_ion, n_ion + n_electron + n_beam, q_electron, n_x, dx, c);
        //半整数格子点上に再定義
        for (int i = 0; i < n_x; i++) {
            if (i == n_x-1) {
                current[0][i] = (current_tmp[0][i] + current_tmp[0][0]) / 2.0;
                current[1][i] = current_tmp[1][i];
                current[2][i] = current_tmp[2][i];
            } else {
            current[0][i] = (current_tmp[0][i] + current_tmp[0][i+1]) / 2.0;
            current[1][i] = current_tmp[1][i];
            current[2][i] = current_tmp[2][i];
            }
        }
        //STEP6
        time_evolution_B(B, E, n_x, dx, dt/2.0);
        //STEP7
        time_evolution_E(E, B, current, n_x, dx, dt, c, epsilon0);
        //STEP8
        time_evolution_x(r, gamma, v, 
                         n_particle, dt/2.0, c);
        periodic_boudary_condition_x(r, n_particle, x_max);

    }
}