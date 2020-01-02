import numpy as np
from scipy.interpolate import interp1d
import util

def get_Cgt_limber(lmax, b, h0, omega_m, r, Theta, Dr, Hr, fr, kh, pk):
    """Calculates the cross - angular power spectra for galaxies and ISW using
    the limber small angle approximation.

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
    c = 3e8
    l = []
    Cgt = []
    for i in range(2, lmax+1):
        l.append(float(i))
        interp_p = interp1d(kh, pk)
        P = np.zeros(len(r))
        c1 = np.where(r != 0.)[0]
        condition = c1[np.where((((float(i) + 0.5 )/r[c1]) > kh.min()) & (((float(i) + 0.5 )/r[c1]) < kh.max()))[0]]
        P[condition] = interp_p((float(i) + 0.5)/r[condition])*((Mpc/h0)**3.)
        integrand = Theta*(Dr**2.)*(Hr*1e3/Mpc)*(fr-1.)*P
        Cgt_val = util.integrate(r*Mpc/h0, integrand, total=True)
        Cgt_val *= - 3. * b * ((100.*h0*1e3/Mpc)**2.)* omega_m
        Cgt_val /= (c**3.)*((float(i) + 0.5)**2.)
        Cgt.append(Cgt_val)
    l, Cgt = np.array(l), np.array(Cgt)
    return l, Cgt


def get_Cgg_limber(lmax, b, h0, r, Theta, Dr, kh, pk):
    """Calculates the angular power spectra for galaxies using the limber small
    angle approximation.

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
    l = []
    Cgg = []
    for i in range(2, lmax+1):
        l.append(float(i))
        interp_p = interp1d(kh, pk)
        P = np.zeros(len(r))
        c1 = np.where(r != 0.)[0]
        condition = c1[np.where((((float(i) + 0.5 )/r[c1]) > kh.min()) & (((float(i) + 0.5 )/r[c1]) < kh.max()))[0]]
        P[condition] = interp_p((float(i) + 0.5 )/r[condition])*((Mpc/h0)**3.)
        integrand = ((Theta/(r*Mpc/h0))**2.)*(Dr**2.)*P
        integrand[0] = 0.
        Cgg_val = util.integrate(r*Mpc/h0, integrand, total=True)
        Cgg_val *= b**2.
        Cgg.append(Cgg_val)
    l, Cgg = np.array(l), np.array(Cgg)
    return l, Cgg
