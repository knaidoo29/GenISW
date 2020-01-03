// c++ libraries
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
// GenISW Libraries
#include "genisw.h"

using namespace std;

// Main --------------------------------------------------------------------- //

int main(int argc, char** argv){

  progress prog;
  prog.start("Spherical Bessel Almn to ISW lm");

  // READING PARAMFILE ------------------------------------------------------ //

  prog.start_process("Reading Parameter File", 1);

  const double Mpc = 3.0856e22, c = 3e8;

  string paramfile, whichparam;

  paramfile = argv[1];

  cout << "Paramfile = " << paramfile << endl;

  string Tcmb_str, h0_str, omega_m_str, omega_l_str;
  bool Tcmb_found = false, h0_found = false, omega_m_found = false;
  bool omega_l_found = false;

  // Reading Cosmological Parameters

  whichparam = "Tcmb";
  extract_param(paramfile, whichparam, &Tcmb_str, &Tcmb_found);
  summarise_param(whichparam, Tcmb_str, Tcmb_found);

  whichparam = "h0";
  extract_param(paramfile, whichparam, &h0_str, &h0_found);
  summarise_param(whichparam, h0_str, h0_found);

  whichparam = "omega_m";
  extract_param(paramfile, whichparam, &omega_m_str, &omega_m_found);
  summarise_param(whichparam, omega_m_str, omega_m_found);

  whichparam = "omega_l";
  extract_param(paramfile, whichparam, &omega_l_str, &omega_l_found);
  summarise_param(whichparam, omega_l_str, omega_l_found);

  double Tcmb, h0, omega_m, omega_l;

  Tcmb = stod(Tcmb_str);
  h0 = stod(h0_str);
  omega_m = stod(omega_m_str);
  omega_l = stod(omega_l_str);

  // Reading Ranges

  string Zmin_str, Zmax_str, Rmin_str, Rmax_str, Zmin_isw_str, Zmax_isw_str;
  bool Zmin_found = false, Zmax_found = false, Rmin_found = false, Rmax_found = false;
  bool Zmin_isw_found = false, Zmax_isw_found = false;

  whichparam = "Zmin";
  extract_param(paramfile, whichparam, &Zmin_str, &Zmin_found);
  summarise_param(whichparam, Zmin_str, Zmin_found);

  whichparam = "Zmax";
  extract_param(paramfile, whichparam, &Zmax_str, &Zmax_found);
  summarise_param(whichparam, Zmax_str, Zmax_found);

  whichparam = "Rmin";
  extract_param(paramfile, whichparam, &Rmin_str, &Rmin_found);
  summarise_param(whichparam, Rmin_str, Rmin_found);

  whichparam = "Rmax";
  extract_param(paramfile, whichparam, &Rmax_str, &Rmax_found);
  summarise_param(whichparam, Rmax_str, Rmax_found);

  whichparam = "Z_min_isw";
  extract_param(paramfile, whichparam, &Zmin_isw_str, &Zmin_isw_found);
  summarise_param(whichparam, Zmin_isw_str, Zmin_found);

  whichparam = "Z_max_isw";
  extract_param(paramfile, whichparam, &Zmax_isw_str, &Zmax_isw_found);
  summarise_param(whichparam, Zmax_isw_str, Zmax_isw_found);

  double Zmin, Zmax, Rmin, Rmax, Zmin_isw, Zmax_isw;

  Zmin = stod(Zmin_str);
  Zmax = stod(Zmax_str);
  Rmin = stod(Rmin_str);
  Rmax = stod(Rmax_str);
  Zmin_isw = stod(Zmin_isw_str);
  Zmax_isw = stod(Zmax_isw_str);

  // Resolution

  string Kmax_str, Lmax_str, Nmax_str, Lmax_isw_str, Nmax_isw_str;
  bool Kmax_found = false, Lmax_found = false, Nmax_found = false;
  bool Lmax_isw_found = false, Nmax_isw_found = false;

  whichparam = "Kmax";
  extract_param(paramfile, whichparam, &Kmax_str, &Kmax_found);
  summarise_param(whichparam, Kmax_str, Kmax_found);

  whichparam = "Lmax";
  extract_param(paramfile, whichparam, &Lmax_str, &Lmax_found);
  summarise_param(whichparam, Lmax_str, Lmax_found);

  whichparam = "Nmax";
  extract_param(paramfile, whichparam, &Nmax_str, &Nmax_found);
  summarise_param(whichparam, Nmax_str, Nmax_found);

  whichparam = "L_max_isw";
  extract_param(paramfile, whichparam, &Lmax_isw_str, &Lmax_isw_found);
  summarise_param(whichparam, Lmax_isw_str, Lmax_isw_found);

  whichparam = "N_max_isw";
  extract_param(paramfile, whichparam, &Nmax_isw_str, &Nmax_isw_found);
  summarise_param(whichparam, Nmax_isw_str, Nmax_isw_found);

  double Kmax, Lmax, Nmax, Lmax_isw, Nmax_isw;

  Kmax = stod(Kmax_str);
  Lmax = stoi(Lmax_str);
  Nmax = stoi(Nmax_str);
  Lmax_isw = stoi(Lmax_isw_str);
  Nmax_isw = stoi(Nmax_isw_str);

  // Inputs

  string inpath, map_roots_LM, map_Redges;
  bool inpath_found = false, map_roots_LM_found = false, map_Redges_found = false;

  whichparam = "inpath";
  extract_param(paramfile, whichparam, &inpath, &inpath_found);
  summarise_param(whichparam, inpath, inpath_found);

  whichparam = "map_roots_LM";
  extract_param(paramfile, whichparam, &map_roots_LM, &map_roots_LM_found);
  summarise_param(whichparam, map_roots_LM, map_roots_LM_found);

  whichparam = "map_Redges";
  extract_param(paramfile, whichparam, &map_Redges, &map_Redges_found);
  summarise_param(whichparam, map_Redges, map_Redges_found);

  // Outputs

  string outpath, sbt_coef_fname, isw_lm_fname;
  bool outpath_found = false, sbt_coef_fname_found = false, isw_lm_fname_found = false;

  whichparam = "outpath";
  extract_param(paramfile, whichparam, &outpath, &outpath_found);
  summarise_param(whichparam, outpath, outpath_found);

  whichparam = "sbt_coef_fname";
  extract_param(paramfile, whichparam, &sbt_coef_fname, &sbt_coef_fname_found);
  summarise_param(whichparam, sbt_coef_fname, sbt_coef_fname_found);

  whichparam = "isw_lm_fname";
  extract_param(paramfile, whichparam, &isw_lm_fname, &isw_lm_fname_found);
  summarise_param(whichparam, isw_lm_fname, isw_lm_fname_found);

  // READING DATA ----------------------------------------------------------- //

  prog.start_process("Loading Spherical Bessel Coefficients Almn ", 1);

  cout << "Loading L and M modes.\n" << endl;

  vector<double>LM;
  int row_length, column_length;

  row_length = read_ascii_table(inpath + map_roots_LM, LM);
  column_length = LM.size()/row_length;

  vector<double>L_double;
  vector<double>M_double;

  extract_from_table(LM, 0, row_length, L_double);
  extract_from_table(LM, 1, row_length, M_double);

  vector<int>L;
  vector<int>M;
  for(int i = 0; i < L_double.size(); i++){
    L.push_back((int)L_double[i]);
    M.push_back((int)M_double[i]);
  }

  cout << "Loading delta_lmn.\n" << endl;

  vector<double>load_Almn;
  vector<double>delta_lmn_real;
  vector<double>delta_lmn_imag;

  row_length = read_ascii_table(outpath + sbt_coef_fname, load_Almn);
  column_length = load_Almn.size()/row_length;

  extract_from_table(load_Almn, 4, row_length, delta_lmn_real);
  extract_from_table(load_Almn, 5, row_length, delta_lmn_imag);

  print_vector(delta_lmn_imag);

  // Spherical Bessel Coefficients Almn to ISW Alm -------------------------- //

  prog.start_process("Spherical Bessel Coefficiens Almn to ISW Alm", 1);

  vector<double>isw_lm_real;
  vector<double>isw_lm_imag;
  double qnl, knl, Nnl, Rnl1, Rnl2, f1, f2;
  double isw_lm_val, isw_lm_real_val, isw_lm_imag_val;

  int zNum = 100;
  vector<double>zz;
  vector<double>rr;
  double dz, dr;
  dz = (Zmax_isw-Zmin_isw)/((double)(zNum-1));
  for(int i = 0; i < zNum; i++){
    zz.push_back(Zmin_isw+dz*((double)i));
  }

  for(int i = 0; i < zz.size(); i++){
    rr.push_back(get_r(zz[i], h0, omega_m, omega_l));
  }

  double r1, r2, H1, H2, fz1, fz2, D0, D1, D2;

  D0 = get_D_unnormed(0., h0, omega_m, omega_l);

  vector<double>Dz;

  for(int j = 0; j < zz.size(); j++){
    Dz.push_back(get_D(zz[j], h0, omega_m, omega_l, D0));
  }

  for(int i = 0; i < L.size(); i++){
    if(L[i] < 2){
      isw_lm_real.push_back(0.);
      isw_lm_imag.push_back(0.);
    }
    else{
      isw_lm_real_val = 0.;
      isw_lm_imag_val = 0.;
      for(int n = 1; n <= Nmax; n++){
        qnl = get_qnl(L[i], n);
        knl = get_knl(L[i], n, Rmax*Mpc/h0, qnl);
        Nnl = get_Nnl(L[i], n, Rmax*Mpc/h0, knl);
        isw_lm_val = 0.;
        if(knl <= Kmax*h0/Mpc && L[i] <= Lmax_isw && n <= Nmax_isw){
          for(int j = 1; j < zz.size(); j++){
            r1 = rr[j-1]*Mpc/h0;
            r2 = rr[j]*Mpc/h0;
            dr = r2 - r1;
            D1 = Dz[j-1];
            D2 = Dz[j];
            H1 = get_H(zz[j-1], h0, omega_m, omega_l)*1e3/Mpc;
            H2 = get_H(zz[j], h0, omega_m, omega_l)*1e3/Mpc;
            Rnl1 = get_Rnl(L[i], r1, knl, Nnl);
            Rnl2 = get_Rnl(L[i], r2, knl, Nnl);
            fz1 = get_fz(zz[j-1], omega_m, omega_l);
            fz2 = get_fz(zz[j], omega_m, omega_l);
            f1 = D1*H1*(1.-fz1)*Rnl1;
            f2 = D2*H2*(1.-fz2)*Rnl2;
            isw_lm_val += 0.5*dr*(f1+f2);
          }
        }
        isw_lm_real_val += (delta_lmn_real[(Nmax*i)+n-1]/pow(knl, 2.))*isw_lm_val;
        isw_lm_imag_val += (delta_lmn_imag[(Nmax*i)+n-1]/pow(knl, 2.))*isw_lm_val;
      }
      isw_lm_real_val *= 3.*omega_m*(pow(100.*h0*1e3/Mpc, 2.)/pow(c, 3.));
      isw_lm_imag_val *= 3.*omega_m*(pow(100.*h0*1e3/Mpc, 2.)/pow(c, 3.));
      isw_lm_real.push_back(isw_lm_real_val*Tcmb);
      isw_lm_imag.push_back(isw_lm_imag_val*Tcmb);
    }
    cout << "delta_lmn to ISW_lm : "<< i+1 << " / " << L.size() << endl;
  }

  // Save Data -------------------------------------------------------------- //

  prog.start_process("Saving ISW lm", 1);

  cout << "Saving to file : " << outpath + isw_lm_fname << endl;

  writedat(outpath + isw_lm_fname, L.begin(), L.end(), M.begin(), M.end(),
    isw_lm_real.begin(), isw_lm_real.end(), isw_lm_imag.begin(), isw_lm_imag.end());

  prog.end();

  return 0;
}
