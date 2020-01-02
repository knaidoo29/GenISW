# GenISW

**Currently being tested and therefore should not be used yet.**

The main goal for GenISW is to generate Integrated Sachs-Wolfe maps from healpix density contrast maps provided in radial bins.

The program is broken into two parts:

1. GenISW (C++): C++ routine for generating Integrate Sachs-Wolfe using Spherical Bessel transform.
  * SH2SBT: Converts a the spherical harmonics for a set of maps into their Spherical Bessel coefficients.
  * SBT2ISW: Outputs the ISW spherical harmonics from the spherical bessel coefficients of the density contrast.
2. pyGenISW (Python): Python function.
  * Calculates the theoretical auto and cross angular power spectra for the density field and ISW in the linear regime
  * Generates parameter files and input files for the C++ scripts.

## Dependencies

GenISW:

* OpenMPI
* GSL
* boost

pyGenISW:

* numpy
* scipy
* camb
* healpy

## Installation

Edit the makefile in the main GenISW directory to install GenISW library which contains the main library of functions. Then edit the makefile in the scripts directory to install the SH2SBT, SH2SBT_MPI, SBT2ISW and SBT2ISW_MPI scripts. pyGenISW can be installed by simply adding the path to GenISW to your pythonpath.
