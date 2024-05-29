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
#include <limits>

using namespace std;


///////////////////////////////////////////////////////////////////////


int main()
{
    string dirname = "mr2008_d=0.01", filename = "mr2008";
    //if(!filesystem::exists(dirname)) {
    //    filesystem::create_directory(dirname);
    //}
    double inf = numeric_limits<double>::infinity();
    int inf_int = numeric_limits<int>::infinity();

    const double c = 0.5, epsilon0 = 1.0, mu_0 = 1.0 / (epsilon0 * pow(c, 2));
    int n_e = 10, n_i = 10, step = 3000;
    int n_x, n_y;
    int n_electron, n_electron_background, n_ion, n_ion_background, n_particle;
    const int iteration = 10000; 
    double dx = 1.0, dy = 1.0, x_max = n_x * dx, y_max = n_y * dy, dt = 1.0;
    int buffer_ion, buffer_electron;
    int total_ion, total_electron;
    int in_ion = 0, in_electron = 0, out_ion = 0, out_electron = 0;
    int in_ion1 = 0, in_electron1 = 0, out_ion1 = 0, out_electron1 = 0;
    int in_ion2 = 0, in_electron2 = 0, out_ion2 = 0, out_electron2 = 0;
    double dmin = 0.01, dmax = 0.01; //補正用。A.B.Langdon(1992)

    double B0, B0_g;
    double t_r, T_e, T_i;
    double m_unit, r_m, m_electron, m_ion;
    double q_unit, r_q, q_ion, q_electron;
    double omega_pe, omega_pi, omega_ce, omega_ci;
    double V_A, C_S, debye_length;
    double ion_inertial_length, sheat_thickness;
    vector<double> v_ion(3, 0), v_electron(3, 0);
    double v_thi, v_the, vb_thi, vb_the;

    m_unit = 1.0;
    r_m = 1.0/4.0;
    t_r = 1.0;
    m_electron = 1 * m_unit;
    m_ion = m_electron / r_m;
    r_q = 1.0;
    B0 = sqrt(static_cast<double>(n_e)) / 1.0;
    T_i = (B0 * B0 / 2.0 / mu_0) / (n_i + n_e * t_r);
    T_e = T_i * t_r;
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
    ion_inertial_length = c / omega_pi;
    sheat_thickness = 2.0 * ion_inertial_length;

    dt = 1.0;
    dx = debye_length;
    dy = debye_length;
    n_x = int(ion_inertial_length * 70);
    n_y = int(ion_inertial_length * 35);
    x_max = n_x * dx - dx;
    y_max = n_y * dy - dy;
    v_thi = sqrt(T_i / m_ion);
    v_the = sqrt(T_e / m_electron);
    vb_thi = sqrt(T_i / 10.0 / m_ion);
    vb_the = sqrt(T_e / 10.0 / m_electron);
    v_electron[0] = 0.0;
    v_electron[1] = 0.0;
    v_electron[2] = c * debye_length / sheat_thickness * sqrt(2 / (1.0 + 1.0/t_r));
    v_ion[0] = -v_electron[0] / t_r;
    v_ion[1] = -v_electron[1] / t_r;
    v_ion[2] = -v_electron[2] / t_r;

    n_ion = int(n_x * n_i * 2.0 * sheat_thickness);
    n_electron = int(n_ion * r_q);
    n_ion_background = int(n_x * 0.2 * n_i * (y_max - 2.0 * sheat_thickness));
    n_electron_background = int(n_x * 0.2 * n_e * (y_max - 2.0 * sheat_thickness));
    buffer_ion = n_i * 10 * n_y; //outflowだし10グリッド分くらい取っておけば大丈夫だろう
    buffer_electron = n_e * 10 * n_y; //outflowだし10グリッド分くらい取っておけば大丈夫だろう
    //シートイオン・背景イオン・バッファーイオン・シート電子・背景電子・バッファー電子の形の配列にする。
    n_particle = n_ion + n_ion_background + buffer_ion + n_electron + n_electron_background + buffer_electron;
    total_ion = n_ion + n_ion_background; //バッファーは入れていない。コード内で変化するパラメータ。
    total_electron = n_electron + n_electron_background; //バッファーは入れていない。コード内で変化するパラメータ。

    cout << "Total number of particles is " << n_particle << ".\n";
    
///////////////////////////////////////////////////////////////////////

    vector<vector<vector<double>>> B(3, vector<vector<double>>(n_x, vector<double>(n_y, 0.0)));
    vector<vector<vector<double>>> B_tmp(3, vector<vector<double>>(n_x, vector<double>(n_y, 0.0)));
    vector<vector<vector<double>>> E(3, vector<vector<double>>(n_x, vector<double>(n_y, 0.0)));
    vector<vector<vector<double>>> E_tmp(3, vector<vector<double>>(n_x, vector<double>(n_y, 0.0)));
    vector<vector<vector<double>>> current(3, vector<vector<double>>(n_x, vector<double>(n_y, 0.0)));
    vector<vector<vector<double>>> current_tmp(3, vector<vector<double>>(n_x, vector<double>(n_y, 0.0)));
    vector<vector<double>> rho(n_x, vector<double>(n_y, 0.0));
    vector<vector<double>> delta_rho(n_x, vector<double>(n_y, 0.0));
    vector<vector<double>> F(n_x, vector<double>(n_y, 0.0));
    vector<vector<double>> phi(n_x, vector<double>(n_y, 0.0));
    vector<vector<double>> zeroth_moment_ion(n_x, vector<double>(n_y, 0.0));
    vector<vector<double>> zeroth_moment_electron(n_x, vector<double>(n_y, 0.0));
    vector<vector<vector<double>>> first_moment_ion(3, vector<vector<double>>(n_x, vector<double>(n_y, 0.0)));
    vector<vector<vector<double>>> first_moment_electron(3, vector<vector<double>>(n_x, vector<double>(n_y, 0.0)));
    vector<vector<vector<double>>> second_moment_ion(9, vector<vector<double>>(n_x, vector<double>(n_y, 0.0)));
    vector<vector<vector<double>>> second_moment_electron(9, vector<vector<double>>(n_x, vector<double>(n_y, 0.0)));
    //バグとりのため、infでセットする
    vector<double> r(3 * n_particle, inf);
    vector<double> r_past(3 * n_particle, inf);
    vector<double> v(3 * n_particle, inf);
    vector<double> v_past(3 * n_particle, inf);
    vector<double> gamma(n_particle, 0.0);
    vector<double> B_particle(3 * n_particle, 0.0);
    vector<double> E_particle(3 * n_particle, 0.0);
    vector<int> index_ion(n_particle/2, 0);
    vector<int> index_electron(n_particle/2, 0);
    vector<int> r_index(3 * n_particle, 0);
    vector<double> cr(3 * n_particle, 0.0);
    vector<double> v_minus(3 * n_particle, 0.0);
    vector<double> v_0(3 * n_particle, 0.0);
    vector<double> v_plus(3 * n_particle, 0.0);
    vector<double> S(3 * n_particle, 0.0);
    vector<double> T(3 * n_particle, 0.0);
    vector<double> tmp_copy_double(n_particle, 0.0);
    vector<double> in_r_ion(3 * buffer_ion, 0.0);
    vector<double> in_v_ion(3 * buffer_ion, 0.0);
    vector<double> in_r_electron(3 * buffer_electron, 0.0);
    vector<double> in_v_electron(3 * buffer_electron, 0.0);
    vector<double> out_r_ion(3 * buffer_ion, 0.0);
    vector<double> out_v_ion(3 * buffer_ion, 0.0);
    vector<double> out_r_electron(3 * buffer_electron, 0.0);
    vector<double> out_v_electron(3 * buffer_electron, 0.0);
    vector<double> cross_gamma_ion(buffer_ion, 0.0);
    vector<double> cross_gamma_electron(buffer_electron, 0.0);
    vector<int> cross_r_index(3 * buffer_ion, 0); // buffer_electronでもいい。
    vector<double> cross_cr(3 * buffer_ion, 0.0); // buffer_electronでもいい。

    ///////////////////////////////////////////////////////////////////////

    //反平行磁場とガイド磁場の設定
    B0_g = 0.0 * B0;
    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_y; j++) {
            B[0][i][j] = B0 * tanh((j*dy - y_max/2.0) / sheat_thickness);
            B[2][i][j] = B0_g;
        }
    }

    //トリガーを先に加えておくことにする
    const double reconnection_ratio = 0.1, Xpoint_position = 10.0 * ion_inertial_length;
    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_y; j++) {
            B[0][i][j] += -B0 * reconnection_ratio * (j*dy - y_max/2.0) / sheat_thickness
                        * exp(-(pow((i*dx - Xpoint_position), 2) + pow((j*dy - y_max/2.0), 2))
                        / pow(2.0 * sheat_thickness, 2));
            B[1][i][j] += B0 * reconnection_ratio * (i*dx - Xpoint_position) / sheat_thickness
                        * exp(-(pow((i*dx - Xpoint_position), 2) + pow((j*dy - y_max/2.0), 2))
                        / pow(2.0 * sheat_thickness, 2));   
        }
    }

    //メンドウなので、n_ion+n_ion_background+n_electron+n_electron_backgroundはn_particleで
    set_initial_position_x(r, 0, total_ion, 0.5*dx, x_max, 1);
    set_initial_position_y_harris(r, 0, n_ion, 0.5*dy, y_max, sheat_thickness, 2);
    set_initial_position_y_background(r, n_ion, total_ion, 0.5*dy, y_max, sheat_thickness, 3);
    set_initial_position_x(r, n_particle/2, n_particle/2+total_electron, 0.5*dx, x_max, 4);
    set_initial_position_y_harris(r, n_particle/2, n_particle/2+n_electron, 0.5*dy, y_max, sheat_thickness, 5);
    set_initial_position_y_background(r, n_particle/2+n_electron, n_particle/2+total_electron, 0.5*dy, y_max, sheat_thickness, 6);
    set_initial_position_z(r, 0, total_ion);
    set_initial_position_z(r, n_particle/2, n_particle/2+total_electron);
    set_initial_velocity_x(v, 0, n_ion, v_ion[0], v_thi, 100); 
    set_initial_velocity_y(v, 0, n_ion, v_ion[1], v_thi, 101); 
    set_initial_velocity_z(v, 0, n_ion, v_ion[2], v_thi, 102); 
    set_initial_velocity_x(v, n_ion, total_ion, 0.0, vb_thi, 103); 
    set_initial_velocity_y(v, n_ion, total_ion, 0.0, vb_thi, 104); 
    set_initial_velocity_z(v, n_ion, total_ion, 0.0, vb_thi, 105); 
    set_initial_velocity_x(v, n_particle/2, n_particle/2+n_electron, v_electron[0], v_the, 106); 
    set_initial_velocity_y(v, n_particle/2, n_particle/2+n_electron, v_electron[1], v_the, 107); 
    set_initial_velocity_z(v, n_particle/2, n_particle/2+n_electron, v_electron[2], v_the, 108); 
    set_initial_velocity_x(v, n_particle/2+n_electron, n_particle/2+total_electron, 0.0, vb_the, 109); 
    set_initial_velocity_y(v, n_particle/2+n_electron, n_particle/2+total_electron, 0.0, vb_the, 110); 
    set_initial_velocity_z(v, n_particle/2+n_electron, n_particle/2+total_electron, 0.0, vb_the, 111);  

    sort_particles(index_ion, r, v, 
                   tmp_copy_double, 0, total_ion, n_x, n_y);
    sort_particles(index_electron, r, v, 
                   tmp_copy_double, n_particle/2, n_particle/2+total_electron, n_x, n_y);

    get_rho(rho, r_index, cr, r, 
            0, total_ion, q_ion, n_x, n_y, dx, dy);
    get_rho(rho, r_index, cr, r, 
            n_particle/2, n_particle/2+total_electron, q_electron, n_x, n_y, dx, dy);
    get_rho_open(rho, cross_r_index, cross_cr, out_r_ion, 
                    0, in_ion, q_ion, n_x, n_y, dx, dy);
    get_rho_open(rho, cross_r_index, cross_cr, out_r_electron, 
                    0, in_electron, q_electron, n_x, n_y, dx, dy);
    poisson_solver_jacobi(phi, rho, iteration, n_x, n_y, dx, dy, epsilon0);
    get_E(E, phi, n_x, n_y, dx, dy);
    boundary_E(E, rho, n_x, n_y, dx, dy, epsilon0);
    filter_E(E, rho, F, n_x, n_y, dx, dy, dt, epsilon0, dmin, dmax);

    double tmp4 = 1.0 / (c * c);
    for (int i = 0; i < 3 * total_ion; i+=3) {
        gamma[i/3] = sqrt(1.0 + (v[i] * v[i] + v[i+1] * v[i+1] + v[i+2] * v[i+2]) * tmp4);
    }
    for (int i = 3 * n_particle / 2; i < 3 * (n_particle / 2 + total_electron); i+=3) {
        gamma[i/3] = sqrt(1.0 + (v[i] * v[i] + v[i+1] * v[i+1] + v[i+2] * v[i+2]) * tmp4);
    }

    string record;
    record = "step.txt";
    ofstream file_record(record);
    file_record.open(record, ios::app);
    file_record << int(0) << "step done..." << "\n";
    file_record.close();


    for (int k = 0; k < step+1; k++) {

        //if (k == int(100/omega_ci)) {
        //    dmin *= 5; 
        //    dmax *= 5;
        //}

        if (k % 10 == 0) {
            file_record.open(record, ios::app);
            file_record << int(k*dt) << "step done..." << "\n";
            file_record.close();
            cout << setw(5) << int(k*dt) << "step done...  " 
                 << " ion : " << setw(8) << total_ion 
                 << " electron : " << setw(8) <<  total_electron
                 << " in ion : " << setw(3) << in_ion1 + in_ion2 
                 << " in electron : " << setw(3) << in_electron1 + in_electron2 
                 << " out ion : " << setw(3) << out_ion1 + out_ion2 
                 << " out electron : " << setw(3) << out_electron1 + out_electron2
                 << " d = " << setw(3) << dmin << "\n";
        }

        //OUTPUT
        if (k % 100 == 0) {
            io_xvKE(r, v, 
                    out_r_ion, out_v_ion, 
                    out_r_electron, out_v_electron, 
                    dirname, filename, k, 
                    total_ion, total_electron, n_particle, 
                    in_ion1 + in_ion2, in_electron1 + in_electron2, 
                    out_ion1 + out_ion2, out_electron1 + out_electron2, 
                    n_x, n_y, m_ion, m_electron);
        }   
        if (k % 100 == 0) {
            for (int i = 0; i < n_x; i++) {
                for (int j = 0; j < n_y; j++) {
                    zeroth_moment_ion[i][j] = 0.0;
                    zeroth_moment_electron[i][j] = 0.0;
                }
            }
            for (int comp = 0; comp < 3; comp++) {
                for (int i = 0; i < n_x; i++) {
                    for (int j = 0; j < n_y; j++) {
                        first_moment_ion[comp][i][j] = 0.0;
                        first_moment_electron[comp][i][j] = 0.0;
                    }
                }
            }
            for (int comp = 0; comp < 9; comp++) {
                for (int i = 0; i < n_x; i++) {
                    for (int j = 0; j < n_y; j++) {
                        second_moment_ion[comp][i][j] = 0.0;
                        second_moment_electron[comp][i][j] = 0.0;
                    }
                }
            }
            get_moment(zeroth_moment_ion, first_moment_ion, second_moment_ion, 
                       r_index, cr, gamma, r, v, 0, total_ion, n_x, n_y, dx, dy, c);
            get_moment(zeroth_moment_electron, first_moment_electron, second_moment_electron, 
                       r_index, cr, gamma, r, v, n_particle/2, n_particle/2+total_electron, n_x, n_y, dx, dy, c);
            io_EBmoment(E, B, rho, 
                        zeroth_moment_ion, zeroth_moment_electron, 
                        first_moment_ion, first_moment_electron, 
                        second_moment_ion, second_moment_electron, 
                        dirname, filename, k, 
                        n_x, n_y, dx, dy, epsilon0, mu_0);
        }     
        //STEP2
        time_evolution_B(B, E, n_x, n_y, dx, dy, dt/2.0);
        boundary_B(B, n_x, n_y, dx, dy, B0_g);

        //STEP3
        //電場・磁場を整数格子点上に再定義
        for (int i = 1; i < n_x; i++) {
            for (int j = 1; j < n_y; j++) {
                B_tmp[0][i][j] = (B[0][i][j] + B[0][i][(n_y+j-1)%n_y]) / 2.0;
                B_tmp[1][i][j] = (B[1][i][j] + B[1][(n_x+i-1)%n_x][j]) / 2.0;
                B_tmp[2][i][j] = (B[2][i][j] + B[2][(n_x+i-1)%n_x][j] + B[2][i][(n_y+j-1)%n_y] + B[2][(n_x+i-1)%n_x][(n_y+j-1)%n_y]) / 4.0;
                E_tmp[0][i][j] = (E[0][i][j] + E[0][(n_x+i-1)%n_x][j]) / 2.0;
                E_tmp[1][i][j] = (E[1][i][j] + E[1][i][(n_y+j-1)%n_y]) / 2.0;
                E_tmp[2][i][j] = E[2][i][j];
            }
        }
        for (int i = 0; i < n_x; i++) {
            B_tmp[0][i][0] = B_tmp[0][i][1];
            B_tmp[1][i][0] = 0.0;
            B_tmp[2][i][0] = B_tmp[2][i][1];
        }
        for (int i = 0; i < n_x; i++) {
            E_tmp[0][i][0] = 0.0;
            E_tmp[1][i][0] = E_tmp[1][i][1];
            E_tmp[2][i][0] = 0.0;
        }
        for (int j = 0; j < n_y; j++) {
            B_tmp[0][0][j] = B_tmp[0][1][j];
            B_tmp[1][0][j] = 0.0;
            B_tmp[2][0][j] = B0_g;
        }
        for (int j = 0; j < n_y; j++) {
            E_tmp[0][0][j] = 0.0;
            E_tmp[1][0][j] = E_tmp[1][1][j];
            E_tmp[2][0][j] = E_tmp[2][1][j];
        }
        fill(B_particle.begin(), B_particle.end(), 0.0);
        fill(E_particle.begin(), E_particle.end(), 0.0);
        get_particle_field(B_particle, E_particle, 
                           r_index, cr,  
                           B_tmp, E_tmp, r, 0, total_ion, n_x, n_y, dx, dy);
        get_particle_field(B_particle, E_particle, 
                           r_index, cr, 
                           B_tmp, E_tmp, r, n_particle/2, n_particle/2+total_electron, n_x, n_y, dx, dy);
        copy(v.begin(), v.end(), v_past.begin()); 
        time_evolution_v(v, gamma, T, S, v_minus, v_0, v_plus, 
                         B_particle, E_particle, 
                         0, total_ion, 
                         m_ion, q_ion, dt, c); 
        time_evolution_v(v, gamma, T, S, v_minus, v_0, v_plus, 
                         B_particle, E_particle, 
                         n_particle/2, n_particle/2+total_electron, 
                         m_electron, q_electron, dt, c); 

        //STEP4
        copy(r.begin(), r.end(), r_past.begin());
        time_evolution_x(r, gamma, v, 0, total_ion, dt/2.0, c);
        time_evolution_x(r, gamma, v, n_particle/2, n_particle/2+total_electron, dt/2.0, c);
        refrective_boudary_condition_x_left(v, r, 0, total_ion, dx/2.0);
        refrective_boudary_condition_x_left(v, r, n_particle/2, n_particle/2+total_electron, dx/2.0);
        open_boudary_condition_x_right(v, r, r_past, in_r_ion, in_r_electron, out_r_ion, out_v_ion, in_ion, out_ion, 
                                       0, total_ion, n_particle, x_max, dx);
        open_boudary_condition_x_right(v, r, r_past, in_r_electron, in_v_electron, out_r_electron, out_v_electron, in_electron, out_electron, 
                                       n_particle/2, n_particle/2+total_electron, n_particle, x_max, dx);      
        
        sort_buffer(v, r, v_past, r_past, total_ion, total_electron, n_particle);
        tmp4 = 1.0 / (c * c);
        for (int i = 0; i < 3 * total_ion; i+=3) {
            gamma[i/3] = sqrt(1.0 + (v[i] * v[i] + v[i+1] * v[i+1] + v[i+2] * v[i+2]) * tmp4);
        }
        for (int i = 3 * n_particle / 2; i < 3 * (n_particle / 2 + total_electron); i+=3) {
            gamma[i/3] = sqrt(1.0 + (v[i] * v[i] + v[i+1] * v[i+1] + v[i+2] * v[i+2]) * tmp4);
        }
        for (int i = 0; i < 3 * out_ion; i+=3) {
            cross_gamma_ion[i/3] = sqrt(1.0 + (out_v_ion[i] * out_v_ion[i] + out_v_ion[i+1] * out_v_ion[i+1] + out_v_ion[i+2] * out_v_ion[i+2]) * tmp4);
        }
        for (int i = 0; i < 3 * out_electron; i+=3) {
            cross_gamma_electron[i/3] = sqrt(1.0 + (out_v_electron[i] * out_v_electron[i] + out_v_electron[i+1] * out_v_electron[i+1] + out_v_electron[i+2] * out_v_electron[i+2]) * tmp4);
        }
        refrective_boudary_condition_y(v, r, 0, total_ion, dy/2.0, y_max);
        refrective_boudary_condition_y(v, r, n_particle/2, n_particle/2+total_electron, dy/2.0, y_max);
        in_ion1 = in_ion; in_electron1 = in_electron; out_ion1 = out_ion; out_electron1 = out_electron;

        //STEP5
        //0で初期化
        for (int comp = 0; comp < 3; comp++) {
            for (int i = 0; i < n_x; i++) {
                for (int j = 0; j < n_y; j++) {
                    current_tmp[comp][i][j] = 0.0;
                }
            }
        }
        get_current_density(current_tmp, r_index, 
                            cr, gamma, r, v, 
                            0, total_ion, q_ion, n_x, n_y, dx, dy, c);
        get_current_density(current_tmp, r_index, 
                            cr, gamma, r, v, 
                            n_particle/2, n_particle/2+total_electron, q_electron, n_x, n_y, dx, dy, c);
        get_current_density_open(current_tmp, cross_r_index, 
                                 cross_cr, cross_gamma_ion, out_r_ion, out_v_ion, 
                                 0, out_ion, q_ion, n_x, n_y, dx, dy, c);
        get_current_density_open(current_tmp, cross_r_index, 
                                 cross_cr, cross_gamma_electron, out_r_electron, out_v_electron, 
                                 0, out_electron, q_electron, n_x, n_y, dx, dy, c);
        //半整数格子点上に再定義
        for (int i = 0; i < n_x-1; i++) {
            for (int j = 0; j < n_y; j++) {
                current[0][i][j] = (current_tmp[0][i][j] + current_tmp[0][(i+1)%n_x][j]) / 2.0;
                current[1][i][j] = (current_tmp[1][i][j] + current_tmp[1][i][(j+1)%n_y]) / 2.0;
                current[2][i][j] = current_tmp[2][i][j];
            }
        }
        //とりあえずn_x-1で代用
        for (int j = 0; j < n_y; j++) {
            current[0][n_x-1][j] = current_tmp[0][n_x-2][j];
            current[1][n_x-1][j] = current_tmp[1][n_x-1][j];
            current[2][n_x-1][j] = current_tmp[2][n_x-1][j];
        }

        //STEP6
        time_evolution_B(B, E, n_x, n_y, dx, dy, dt/2.0);
        boundary_B(B, n_x, n_y, dx, dy, B0_g);
        
        //STEP7
        for (int i = 0; i < n_x; i++) {
            for (int j = 0; j < n_y; j++) {
                rho[i][j] = 0.0;
            }
        }
        get_rho(rho, r_index, cr, r, 
                0, total_ion, q_ion, n_x, n_y, dx, dy);
        get_rho(rho, r_index, cr, r, 
                n_particle/2, n_particle/2+total_electron, q_electron, n_x, n_y, dx, dy);
        get_rho_open(rho, cross_r_index, cross_cr, out_r_ion, 
                     0, out_ion, q_ion, n_x, n_y, dx, dy);
        get_rho_open(rho, cross_r_index, cross_cr, out_r_electron, 
                     0, out_electron, q_electron, n_x, n_y, dx, dy);
        time_evolution_E(E, B, current, n_x, n_y, dx, dy, dt, c, epsilon0);
        boundary_E(E, rho, n_x, n_y, dx, dy, epsilon0);
        filter_E(E, rho, F, n_x, n_y, dx, dy, dt, epsilon0, dmin, dmax);
        
        //STEP8
        copy(r.begin(), r.end(), r_past.begin());
        time_evolution_x(r, gamma, v, 0, total_ion, dt/2.0, c);
        time_evolution_x(r, gamma, v, n_particle/2, n_particle/2+total_electron, dt/2.0, c);
        refrective_boudary_condition_x_left(v, r, 0, total_ion, dx/2.0);
        refrective_boudary_condition_x_left(v, r, n_particle/2, n_particle/2+total_electron, dx/2.0);
        open_boudary_condition_x_right(v, r, r_past, in_r_ion, in_r_electron, out_r_ion, out_v_ion, in_ion, out_ion, 
                                       0, total_ion, n_particle, x_max, dx);
        open_boudary_condition_x_right(v, r, r_past, in_r_electron, in_v_electron, out_r_electron, out_v_electron, in_electron, out_electron, 
                                       n_particle/2, n_particle/2+total_electron, n_particle, x_max, dx);     
        sort_buffer(v, r, v_past, r_past, total_ion, total_electron, n_particle);
        tmp4 = 1.0 / (c * c);
        for (int i = 0; i < 3 * total_ion; i+=3) {
            gamma[i/3] = sqrt(1.0 + (v[i] * v[i] + v[i+1] * v[i+1] + v[i+2] * v[i+2]) * tmp4);
        }
        for (int i = 3 * n_particle / 2; i < 3 * (n_particle / 2 + total_electron); i+=3) {
            gamma[i/3] = sqrt(1.0 + (v[i] * v[i] + v[i+1] * v[i+1] + v[i+2] * v[i+2]) * tmp4);
        }
        for (int i = 0; i < 3 * out_ion; i+=3) {
            cross_gamma_ion[i/3] = sqrt(1.0 + (out_v_ion[i] * out_v_ion[i] + out_v_ion[i+1] * out_v_ion[i+1] + out_v_ion[i+2] * out_v_ion[i+2]) * tmp4);
        }
        for (int i = 0; i < 3 * out_electron; i+=3) {
            cross_gamma_electron[i/3] = sqrt(1.0 + (out_v_electron[i] * out_v_electron[i] + out_v_electron[i+1] * out_v_electron[i+1] + out_v_electron[i+2] * out_v_electron[i+2]) * tmp4);
        }
        refrective_boudary_condition_y(v, r, 0, total_ion, dy/2.0, y_max);
        refrective_boudary_condition_y(v, r, n_particle/2, n_particle/2+total_electron, dy/2.0, y_max);
        in_ion2 = in_ion; in_electron2 = in_electron; out_ion2 = out_ion; out_electron2 = out_electron;


        if (k % 50 == 0) {
            sort_particles(index_ion, r, v, 
                   tmp_copy_double, 0, total_ion, n_x, n_y);
            sort_particles(index_electron, r, v, 
                        tmp_copy_double, n_particle/2, n_particle/2+total_electron, n_x, n_y);
        }
    }
}
