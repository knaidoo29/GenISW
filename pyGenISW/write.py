import numpy as np
import subprocess


def write_paramfile(paramfname, Tcmb, h0, omega_m, omega_l, Zmin, Zmax, Rmin, Rmax,
    Kmax, Lmax, Nmax, inpath, map_roots_LM, map_Redges, outpath, sbt_coef_fname,
    isw_lm_fname, Zmin_isw=None, Zmax_isw=None, Lmax_isw=None, Nmax_isw=None):
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
    Zmin_isw : float, optional
        Sets the minimum redshift for the ISW calculation.
    Zmax_isw : float, optional
        Sets the maximum redshift for the ISW calculation.
    Lmax_isw : float, optional
        Maximum L mode used to calculate the ISW Alm.
    Nmax_isw : float, optional
        Maximum N used to calculate the ISW Almn.
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
    if Zmin_isw is not None:
        paramfile.write("Z_min_isw           " + str(Zmin_isw) + " %\n")
    else:
        paramfile.write("Z_min_isw           " + str(Zmin) + " %\n")
    if Zmax_isw is not None:
        paramfile.write("Z_max_isw           " + str(Zmax_isw) + " %\n")
    else:
        paramfile.write("Z_max_isw           " + str(Zmax) + " %\n")
    if Lmax_isw is not None:
        paramfile.write("L_max_isw           " + str(Lmax_isw) + " %\n")
    else:
        paramfile.write("L_max_isw           " + str(Lmax) + " %\n")
    if Nmax_isw is not None:
        paramfile.write("N_max_isw           " + str(Nmax_isw) + " %\n")
    else:
        paramfile.write("N_max_isw           " + str(Nmax) + " %\n")
    paramfile.close()
