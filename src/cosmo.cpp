#include "../include/genisw.h"
#include <vector>
#include <math.h>

using namespace std;

double get_H(double z, double h0, double omega_m, double omega_l){
  /* Returns the Hubble expansion rate at z.

  Parameters
  ----------
  z : double
    Redshift
  h0 : double
    'Litte' hubble constant.
  omega_m : double
    Matter density.
  omega_l : double
    Lambda density
  */
  return pow(pow(100.*h0, 2.)*(omega_m*pow(1.+z, 3.) + omega_l), 0.5);
}

double get_omega_m_z(double z, double omega_m, double omega_l){
  /* Returns the matter density at z.

  Parameters
  ----------
  z : double
    Redshift
  omega_m : double
    Matter density.
  omega_l : double
    Lambda density
  */
  return (omega_m*pow(1.+z, 3.))/(omega_m*pow(1.+z, 3.) + omega_l);
}

double get_fz(double z, double omega_m, double omega_l){
  /* Returns the differential of the linear growth rate f = dlnD/dlna.

  Parameters
  ----------
  z : double
    Redshift
  omega_m : double
    Matter density.
  omega_l : double
    Lambda density
  */
  return pow(get_omega_m_z(z, omega_m, omega_l), 0.6);
}

double get_r(double z, double h0, double omega_m, double omega_l){
  /* Get comoving distance.

  Parameters
  ----------
  z : double
    Redshift.
  h0 : double
    'Little' hubble constant.
  omega_m : double
    Matter density.
  omega_l : double
    Lambda density.
  */
  double dz, f1, f2, r;
  dz = z/(1000.-1.);
  r = 0.;
  for(int i = 1; i < 1000; i++){
    f1 = 1./pow(omega_m*pow(1. + dz*((double)(i-1)), 3.) + omega_l, 0.5);
    f2 = 1./pow(omega_m*pow(1. + dz*((double)(i)), 3.) + omega_l, 0.5);
    r += 0.5*dz*(f1+f2);
  }
  r *= 3000.;
  return r;
}

double get_D_integrand(double z, double h0, double omega_m, double omega_l){
  /* Returns the integral for calculating the linear growth function.

  Parameters
  ----------
  z : double
    Redshift.
  h0 : double
    'Little' hubble constant.
  omega_m : double
    Matter density.
  omega_l : double
    Lambda density.
  */
  double a;
  a = 1./(1.+z);
  return 1./pow(a*get_H(z, h0, omega_m, omega_l), 3.);
}

double get_D_unnormed(double z, double h0, double omega_m, double omega_l){
  /* Returns the unnormalised linear growth function at redshift z.

  Parameters
  ----------
  z : double
    Redshift.
  h0 : double
    'Little' hubble constant.
  omega_m : double
    Matter density.
  omega_l : double
    Lambda density.
  */
  double a, a_min, da, D_unnormed, z1, z2, f1, f2;
  a = 1./(1.+z);
  a_min = 1e-5;
  da = (a-a_min)/(100.-1.);
  D_unnormed = 0.;
  for(int i = 1; i < 100; i++){
    z1 = (1./(a_min + da*((double)i-1.))) - 1.;
    z2 = (1./(a_min + da*((double)i))) - 1.;
    f1 = get_D_integrand(z1, h0, omega_m, omega_l);
    f2 = get_D_integrand(z2, h0, omega_m, omega_l);
    D_unnormed += 0.5*da*(f1+f2);
  }
  return get_H(z, h0, omega_m, omega_l)*D_unnormed;
}

double get_D(double z, double h0, double omega_m, double omega_l, double D0){
  /* Returns the linear growth function.

  Parameters
  ----------
  z : double
    Redshift.
  h0 : double
    'Little' hubble constant.
  omega_m : double
    Matter density.
  omega_l : double
    Lambda density.
  D0 : double
    The unnormalised linear growth function at z=0.
  */
  return get_D_unnormed(z, h0, omega_m, omega_l)/D0;
}
