//#ifndef COSMO_H
//#define COSMO_H

#include <vector>

using namespace std;

double get_H(double z, double h0, double omega_m, double omega_l);

double get_omega_m_z(double z, double omega_m, double omega_l);

double get_fz(double z, double omega_m, double omega_l);

double get_r(double z, double h0, double omega_m, double omega_l);

double get_D_integrand(double z, double h0, double omega_m, double omega_l);

double get_D_unnormed(double z, double h0, double omega_m, double omega_l);

double get_D(double z, double h0, double omega_m, double omega_l, double D0);

//#endif
