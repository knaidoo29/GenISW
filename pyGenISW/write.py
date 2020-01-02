import numpy as np
import subprocess


def write_paramfile(paramfname, Tcmb, h0, omega_m, omega_l, Zmin, Zmax, Rmin, Rmax,
    Kmax, Lmax, Nmax, inpath, map_roots_LM, map_Redges, outpath, sbt_coef_fname,
    isw_lm_fname, Rmin_isw=None, Rmax_isw=None):
    """ Writes a parameter file to be read by GenISW.

    Parameters
    ----------
    paramfname : str
        Parameter filename.
    Tcmb : float
        CMB temperature.
    h0 : float
        'Little' hubble constant.
    omega_m : float
        Matter density.
    omega_l : float
        Lambda density.
    Zmin : float
        Minimum redshift.
    Zmax : float
        Maximum redshift.
    Rmin : float
        Minimum comoving distance.
    Rmax : float
        Maximum comoving distance.
    Kmax : float
        Maximum k mode.
    Lmax : int
        Maximum l mode.
    Nmax : int
        Maximum n mode.
    inpath : str
        Path to input files.
    map_roots_LM : str
        File defining lm modes.
    map_Redges : str
        Rmin and Rmax for each map.
    outpath : str
        Path to output files.
    sbt_coef_fname : str
        Filename containing SBT coefficients.
    isw_lm_fname : str
        Filename containing the ISW spherical harmonic coefficients.
    Rmin_isw : float, optional
        Sets the minimum distance for the ISW calculation.
    Rmax_isw : float, optional
        Sets the maximum distance for the ISW calculation.
    """
    subprocess.call('touch '+ inpath + paramfname,shell=True)
    paramfile = open(inpath + paramfname, 'w')
    paramfile.write("## Cosmological Parameters ##\n")
    paramfile.write("Tcmb               " + str(Tcmb) + " %\n")
    paramfile.write("h0                 " + str(h0) + " %\n")
    paramfile.write("omega_m            " + str(omega_m) + " %\n")
    paramfile.write("omega_l            " + str(omega_l) + " %\n")
    paramfile.write("## Ranges ##\n")
    paramfile.write("Zmin               " + str(Zmin) + " %\n")
    paramfile.write("Zmax               " + str(Zmax) + " %\n")
    paramfile.write("Rmin               " + str(Rmin) + " %\n")
    paramfile.write("Rmax               " + str(Rmax) + " %\n")
    paramfile.write("## Resolution ##\n")
    paramfile.write("Kmax               " + str(Kmax) + " %\n")
    paramfile.write("Lmax               " + str(Lmax) + " %\n")
    paramfile.write("Nmax               " + str(Nmax) + " %\n")
    paramfile.write("## Inputs ##\n")
    paramfile.write("inpath             " + inpath + " %\n")
    paramfile.write("map_roots_LM       " + map_roots_LM + " %\n")
    paramfile.write("map_Redges         " + map_Redges + " %\n")
    paramfile.write("## Outputs ##\n")
    paramfile.write("outpath            " + outpath  + " %\n")
    paramfile.write("sbt_coef_fname     " + sbt_coef_fname + " %\n")
    paramfile.write("isw_lm_fname       " + isw_lm_fname + " %\n")
    if Rmin_isw is not None:
        paramfile.write("Rmin_isw           " + str(Rmin_isw) + " %\n")
    else:
        paramfile.write("Rmin_isw           " + str(Rmin) + " %\n")
    if Rmin_isw is not None:
        paramfile.write("Rmax_isw           " + str(Rmax_isw) + " %\n")
    else:
        paramfile.write("Rmax_isw           " + str(Rmax) + " %\n")
    paramfile.close()
