#include "../include/genisw.h"
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <boost/math/special_functions/bessel.hpp>

using namespace std;


double get_jl(int l, double r){
  /* Returns the spherical bessel value with mode l at point r using the boost
  libraries.

  Parameters
  ----------
  l : int
    Spherical bessel function mode.
  r : double
    Variable r along the spherical bessel function.
  */
  return boost::math::sph_bessel(l, r);
}


double get_qnl(int l, int n){
  /* Returns the nth root of the spherical bessel function with mode l using the
  GSL root finder for normal bessel functions (i.e. it is shifted by a mode +0.5
  for the spherical bessel root).

  Parameters
  ----------
  l : int
    Spherical bessel function mode.
  n : int
    The nth root.
  */
  return gsl_sf_bessel_zero_Jnu((double)l + 0.5, n);
}

double get_knl(int l, int n, double Rmax, double qnl){
  /* Returns the k mode for the spherical bessel function.

  Parameters
  ----------
  l : int
    Spherical bessel function mode.
  n : int
    The nth root.
  Rmax : double
    The maximum radius.
  qnl : double
    The nth root of the spherical bessel function with mode l.
  */
  return qnl/Rmax;
}

double get_Nnl(int l, int n, double Rmax, double knl){
  /* Returns the spherical bessel transform normalisation constant.

  Parameters
  ----------
  l : int
    Spherical bessel function mode.
  n : int
    The nth root.
  Rmax : double
    The maximum radius.
  knl : double
    The k mode for the spherical bessel function.
  */
  return (pow(Rmax, 3.)/2.)*pow(get_jl(l+1, knl*Rmax), 2.);
}

double get_Rnl(int l, double r, double knl, double Nnl){
  /* The normalised spherical bessel function for spherical bessel transforms.

  Parameters
  ----------
  l : int
    Spherical bessel function mode.
  r : double
    Variable r along the spherical bessel function.
  knl : double
    The k mode for the spherical bessel function.
  Nnl : double
    Spherical bessel transform normalisation constant.
  */
  return get_jl(l, knl*r)/pow(Nnl, 0.5);
}
