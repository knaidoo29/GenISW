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
  prog.start("Spherical Harmonic Alm to Spherical Bessel Coefficients Almn");

  // READING PARAMFILE ------------------------------------------------------ //

  prog.start_process("Reading Parameter File", 1);

  string paramfile, whichparam;

  paramfile = argv[1];

  cout << "Paramfile = " << paramfile << endl;

  string h0_str;
  bool h0_found = false;

  // Reading Cosmological Parameters

  whichparam = "h0";
  extract_param(paramfile, whichparam, &h0_str, &h0_found);
  summarise_param(whichparam, h0_str, h0_found);

  double Mpc, h0, units;

  if(h0_found == true){
    // Assume cosmological units of Mpc/h:
    Mpc = 3.0856e22;
    h0 = stod(h0_str);
    units = Mpc/h0;
  }
  else{
    Mpc = 1.;
    h0 = 1.;
    units = 1.;
  }

  // Reading Ranges

  string Rmin_str, Rmax_str;
  bool Rmin_found = false, Rmax_found = false;

  whichparam = "Rmin";
  extract_param(paramfile, whichparam, &Rmin_str, &Rmin_found);
  summarise_param(whichparam, Rmin_str, Rmin_found);

  whichparam = "Rmax";
  extract_param(paramfile, whichparam, &Rmax_str, &Rmax_found);
  summarise_param(whichparam, Rmax_str, Rmax_found);

  double Rmin, Rmax;

  Rmin = stod(Rmin_str);
  Rmax = stod(Rmax_str);

  // Resolution

  string Kmax_str, Lmax_str, Nmax_str;
  bool Kmax_found = false, Lmax_found = false, Nmax_found = false;

  whichparam = "Kmax";
  extract_param(paramfile, whichparam, &Kmax_str, &Kmax_found);
  summarise_param(whichparam, Kmax_str, Kmax_found);

  whichparam = "Lmax";
  extract_param(paramfile, whichparam, &Lmax_str, &Lmax_found);
  summarise_param(whichparam, Lmax_str, Lmax_found);

  whichparam = "Nmax";
  extract_param(paramfile, whichparam, &Nmax_str, &Nmax_found);
  summarise_param(whichparam, Nmax_str, Nmax_found);

  double Kmax, Lmax, Nmax;

  Kmax = stod(Kmax_str);
  Lmax = stoi(Lmax_str);
  Nmax = stoi(Nmax_str);

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

  string outpath, sbt_coef_fname;
  bool outpath_found = false, sbt_coef_fname_found = false;

  whichparam = "outpath";
  extract_param(paramfile, whichparam, &outpath, &outpath_found);
  summarise_param(whichparam, outpath, outpath_found);

  whichparam = "sbt_coef_fname";
  extract_param(paramfile, whichparam, &sbt_coef_fname, &sbt_coef_fname_found);
  summarise_param(whichparam, sbt_coef_fname, sbt_coef_fname_found);

  // READING DATA ----------------------------------------------------------- //

  prog.start_process("Reading data", 1);

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

  cout << "Loading R edges.\n" << endl;

  vector<double>Redges;

  row_length = read_ascii_table(inpath + map_Redges, Redges);
  column_length = Redges.size()/row_length;

  vector<double>R1;
  vector<double>R2;

  extract_from_table(Redges, 0, row_length, R1);
  extract_from_table(Redges, 1, row_length, R2);

  int NumSlices = R1.size();

  cout << "Loading Alm for slices.\n" << endl;

  vector<double>Alm_slices;
  vector<double>Alm_slices_real;
  vector<double>Alm_slices_imag;

  for(int i = 0; i < NumSlices; i++){

    row_length = read_ascii_table(inpath + "map_alm_" + to_string(i) + ".txt", Alm_slices);
    column_length = Alm_slices.size()/row_length;

    extract_from_table(Alm_slices, 0, row_length, Alm_slices_real);
    extract_from_table(Alm_slices, 1, row_length, Alm_slices_imag);

    Alm_slices.erase(Alm_slices.begin(), Alm_slices.begin()+Alm_slices.size());

    cout << "Reading Map Alm: " + to_string(i) + " / " + to_string(NumSlices) << endl;
  }

  // Spherical Harmonic Alm to Spherical Bessel Coefficients Almn ----------- //

  prog.start_process("Spherical Harmonic Alm to Spherical Bessel Coefficients Almn", 1);

  vector<double>r;

  int rNum = 1000;
  double dr = (Rmax)/((double)rNum-1.);

  for(int i = 0; i < rNum; i++){
    r.push_back(dr*((double)i));
  }

  double qnl, knl, Nnl, Rnl1, Rnl2, f1, f2;

  vector<int>L_mode;
  vector<int>M_mode;
  vector<int>N_mode;
  vector<double>K_mode;
  vector<double>Almn_real;
  vector<double>Almn_imag;
  vector<double>Alm_slice_real;
  vector<double>Alm_slice_imag;

  int index1, index2;
  double Almn_real_val, Almn_imag_val;

  for(int i = 0; i < L.size(); i++){
    for(int n = 1; n <= Nmax; n++){
      qnl = get_qnl(L[i], n);
      knl = get_knl(L[i], n, Rmax*units, qnl);
      Nnl = get_Nnl(L[i], n, Rmax*units, knl);
      L_mode.push_back(L[i]);
      M_mode.push_back(M[i]);
      N_mode.push_back(n);
      K_mode.push_back(knl*units);
      if(knl <= Kmax/units){
        for(int j = 0; j < NumSlices; j++){
          Alm_slice_real.push_back(Alm_slices_real[L.size()*j + i]);
          Alm_slice_imag.push_back(Alm_slices_imag[L.size()*j + i]);
        }
        Almn_real_val = 0.;
        Almn_imag_val = 0.;
        for(int j = 1; j < r.size(); j++){
          Rnl1 = get_Rnl(L[i], r[j-1]*units, knl, Nnl);
          index1 = get_index(r[j-1], R1, R2);
          Rnl2 = get_Rnl(L[i], r[j]*units, knl, Nnl);
          index2 = get_index(r[j], R1, R2);
          f1 = Alm_slice_real[index1]*Rnl1*pow(r[j-1]*units, 2.);
          f2 = Alm_slice_real[index2]*Rnl2*pow(r[j]*units, 2.);
          Almn_real_val += 0.5*(dr*units)*(f1+f2);
          f1 = Alm_slice_imag[index1]*Rnl1*pow(r[j-1]*units, 2.);
          f2 = Alm_slice_imag[index2]*Rnl2*pow(r[j]*units, 2.);
          Almn_imag_val += 0.5*(dr*units)*(f1+f2);
        }
        Almn_real.push_back(Almn_real_val);
        Almn_imag.push_back(Almn_imag_val);
        Alm_slice_real.erase(Alm_slice_real.begin(), Alm_slice_real.begin()+Alm_slice_real.size());
        Alm_slice_imag.erase(Alm_slice_imag.begin(), Alm_slice_imag.begin()+Alm_slice_imag.size());
      }
      else{
        Almn_real.push_back(0.);
        Almn_imag.push_back(0.);
      }
    }
    cout << "Alm to Almn : "<< i+1 << " / " << L.size() << endl;
  }

  // Save Data -------------------------------------------------------------- //

  prog.start_process("Saving Almn", 1);

  cout << "Saving to file : " << outpath + sbt_coef_fname << endl;

  writedat(outpath + sbt_coef_fname,
    L_mode.begin(), L_mode.end(),
    M_mode.begin(), M_mode.end(),
    N_mode.begin(), N_mode.end(),
    K_mode.begin(), K_mode.end(),
    Almn_real.begin(), Almn_real.end(),
    Almn_imag.begin(), Almn_imag.end());

  prog.end();

  return 0;
}
