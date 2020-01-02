//#ifndef SPHERICAL_BESSEL_H
//#define SPHERICAL_BESSEL_H

#include <vector>

using namespace std;

double get_jl(int l, double r);

double get_qnl(int l, int n);

double get_knl(int l, int n, double Rmax, double qnl);

double get_Nnl(int l, int n, double Rmax, double knl);

double get_Rnl(int l, double r, double knl, double Nnl);

//#endif
