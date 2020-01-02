#include "../include/genisw.h"
#include <vector>

using namespace std;

int get_index(double r_val, vector<double> &r1, vector<double> &r2){
  /* Returns index for interpolating within binned values.

  Parameters
  ----------
  r_val : double
    Comoving distance.
  r1 : vector
    Comoving distance on the minimum edge.
  r2 : vector
    Comoving distance on the maximum edge.
  */
  int index;
  for(int i = 0; i < r1.size(); i++){
    if(r_val >= r1[i] and r_val <= r2[i]){
      index = i;
    }
  }
  return index;
}
