import numpy as np
from scipy.interpolate import interp1d
from scipy import special
import util
import angular_power_spec_limber_approx as limber

def pk2delta2k(kh, pk):
    """Convert power spectra to the unitless delta2k.

    Parameters
    ----------
    kh : array
        Power spectrum k given in units of h/Mpc.
    pk : array
        Power spectrum given in units of (Mpc/h)^3.

    Returns
    -------
    delta2k : array
        Unitless delta2k.
    """
    delta2k = (kh**3.) * pk
    delta2k *= 4.*np.pi/((2.*np.pi)**3.)
    return delta2k


def get_Wg(kh, l, r, Theta, Dr, h0):
    """Computes Wg.

    Parameters
    ----------
    kh : array
        Power spectrum k given in units of h/Mpc.
    l : int
        l mode.
    r : array
        Comoving distance given in Mpc/h.
    Theta : array
        Redshift selection function.
    Dr : array
        Linear growth function.
    h0 : float
        'Little' hubble constant.
    """
    integrand = Theta*special.spherical_jn(int(l), kh*r)*Dr
    return util.integrate(r, integrand, total=True)


def get_Wt(kh, l, r, Hr, Dr, fr, h0, omega_m):
    """Computes Wg.

    Parameters
    ----------
    kh : array
        Power spectrum k given in units of h/Mpc.
    l : int
        l mode.
    r : array
        Comoving distance given in Mpc/h.
    Hr : array
        Hubble expansion rate function.
    Dr : array
        Linear growth function.
    fr : array
        Linear growth differential dlnD/dlna.
    h0 : float
        'Little' hubble constant.
    omega_m : float
        Present matter density.
    """
    Mpc = 3.0856e22
    c = 3e8
    integrand = special.spherical_jn(int(l), kh*r)*(Hr)*Dr*(fr - 1.)
    return -1.*3.*omega_m*((100.*h0*1e3/Mpc)**2.)*(1./((kh**2.)*(c**3.)))*util.integrate(r, integrand, total=True)


def get_Cgg_nolimber(lmax, b, h0, r, Theta, Dr, kh, pk):
    """Calculates the angular power spectra for galaxies.

    Parameters
    ----------
    lmax : int
        Maximum l mode.
    b : float
        linear Galaxy bias.
    h0 : float
        'Little' hubble constant.
    r : array
        Comoving distance given in Mpc/h0.
    Theta : array
        Redshift selection function.
    Dr : array
        Linear growth function.
    kh : array
        Power spectrum k given in units of h/Mpc.
    pk : array
        Power spectrum given in units of (Mpc/h)^3.

    Returns
    -------
    l : array
        Angular power spectra l modes.
    Cgg : array
        Angular power spectra for galaxies.
    """
    Mpc = 3.0856e22
    l, Cgg = [], []
    for i in range(2, lmax+1):
        l.append(float(i))
        cgg_integrand = []
        for j in range(0, len(kh)):
            Wg = get_Wg(kh[j]*h0/Mpc, float(i), r*Mpc/h0, Theta, Dr, h0)
            delta2k = pk2delta2k(kh[j]*h0/Mpc, pk[j]*((Mpc/h0)**3))
            cgg_integrand.append(delta2k * (Wg**2.) / (kh[j]*h0/Mpc))
        cgg_integrand = np.array(cgg_integrand)
        Cgg.append(4.*np.pi*(b**2.)*util.integrate(kh*h0/Mpc, cgg_integrand, total=True))
    Cgg = np.array(Cgg)
    l = np.array(l)
    return l, Cgg

def get_Cgg(lmax, b, h0, r, Theta, Dr, kh, pk, lswitch2limber=100):
    """Calculates the angular power spectra for galaxies.

    Parameters
    ----------
    lmax : int
        Maximum l mode.
    b : float
        linear Galaxy bias.
    h0 : float
        'Little' hubble constant.
    r : array
        Comoving distance given in Mpc/h0.
    Theta : array
        Redshift selection function.
    Dr : array
        Linear growth function.
    kh : array
        Power spectrum k given in units of h/Mpc.
    pk : array
        Power spectrum given in units of (Mpc/h)^3.
    lswitch2limber : int, optional
        Switch to limber approximation.

    Returns
    -------
    l : array
        Angular power spectra l modes.
    Cgg : array
        Angular power spectra for galaxies.
    """
    _l, Cgg_full = get_Cgg_nolimber(lswitch2limber, b, h0, r, Theta, Dr, kh, pk)
    l, Cgg_limber = limber.get_Cgg_limber(lmax, b, h0, r, Theta, Dr, kh, pk)
    Cgg = np.copy(Cgg_limber)
    Cgg[:len(Cgg_full)] = Cgg_full
    return l, Cgg


def get_Cgt_nolimber(lmax, b, h0, omega_m, r, Theta, Dr, Hr, fr, kh, pk):
    """Calculates the cross - angular power spectra for galaxies and ISW.

    Parameters
    ----------
    lmax : int
        Maximum l mode.
    b : float
        linear Galaxy bias.
    h0 : float
        'Little' hubble constant.
    omega_m : float
        Present matter density.
    r : array
        Comoving distance given in Mpc/h0.
    Theta : array
        Redshift selection function.
    Dr : array
        Linear growth function.
    Hr : array
        Hubble expansion rate function.
    fr : array
        Linear growth differential dlnD/dlna.
    kh : array
        Power spectrum k given in units of h/Mpc.
    pk : array
        Power spectrum given in units of (Mpc/h)^3.

    Returns
    -------
    l : array
        Angular power spectra l modes.
    Cgt : array
        Angular cross power spectra for galaxies and ISW.
    """
    Mpc = 3.0856e22
    l, Cgt = [], []
    for i in range(2, lmax+1):
        l.append(float(i))
        cgt_integrand = []
        for j in range(0, len(kh)):
            Wg = get_Wg(kh[j]*h0/Mpc, float(i), r*Mpc/h0, Theta, Dr, h0)
            Wt = get_Wt(kh[j]*h0/Mpc, float(i), r*Mpc/h0, Hr*1e3/Mpc, Dr, fr, h0, omega_m)
            delta2k = pk2delta2k(kh[j]*h0/Mpc, pk[j]*(Mpc/h0)**3.)
            cgt_integrand.append(delta2k * (Wg*Wt) / (kh[j]*h0/Mpc))
        cgt_integrand = np.array(cgt_integrand)
        Cgt.append(4.*np.pi*b*util.integrate(kh*h0/Mpc, cgt_integrand, total=True))
    Cgt = np.array(Cgt)
    l = np.array(l)
    return l, Cgt


def get_Cgt(lmax, b, h0, omega_m, r, Theta, Dr, Hr, fr, kh, pk, lswitch2limber=100):
    """Calculates the angular power spectra for galaxies.

    Parameters
    ----------
    lmax : int
        Maximum l mode.
    b : float
        linear Galaxy bias.
    h0 : float
        'Little' hubble constant.
    r : array
        Comoving distance given in Mpc/h0.
    Theta : array
        Redshift selection function.
    Dr : array
        Linear growth function.
    kh : array
        Power spectrum k given in units of h/Mpc.
    pk : array
        Power spectrum given in units of (Mpc/h)^3.
    lswitch2limber : int, optional
        Switch to limber approximation.

    Returns
    -------
    l : array
        Angular power spectra l modes.
    Cgg : array
        Angular power spectra for galaxies.
    """
    _l, Cgt_full = get_Cgt_nolimber(lswitch2limber, b, h0, omega_m, r, Theta, Dr, Hr, fr, kh, pk)
    l, Cgt_limber = limber.get_Cgt_limber(lmax, b, h0, omega_m, r, Theta, Dr, Hr, fr, kh, pk)
    Cgt = np.copy(Cgt_limber)
    Cgt[:len(Cgt_full)] = Cgt_full
    return l, Cgt


def get_Cisw(Cgg, Cgt):
    """Returns the angular power spectra for the ISW:

    Parameters
    ----------
    Cgg : array
        Angular power spectra for galaxies.
    Cgt : array
        Angular cross power spectra for galaxies and ISW.
    """
    return (Cgt**2.)/Cgg
