import numpy as np
from scipy.integrate import quad
import util


def a2z(a):
    """Converts a to z.

    Parameters
    ----------
    a : float
        The scale factor.
    """
    return 1./a - 1.


def z2a(z):
    """Converts z to a.

    Parameters
    ----------
    z : float
        Redshift.
    """
    return 1./(1.+z)


def get_H(a, h0, omega_m, omega_l):
    """Returns the hubble expansion at a.

    Parameters
    ----------
    a : float
        The scale factor.
    h0 : float
        'Little' Hubble constant.
    omega_m : float
        Present matter density.
    omega_l : float
        Present dark energy density
    """
    return np.sqrt(((100.*h0)**2.)*(omega_m*a**-3. + omega_l))


def get_theta_r_constant(r, r_min, r_max, h0):
    """Returns a constant redshift selection function.

    Parameters
    ----------
    r : array
        Comoving distance in  Mpc/h0.
    r_min : float
        Minimum comoving distance.
    r_max : float
        Maximum comoving distance.
    h0 : float
        'Little' Hubble constant.
    """
    Mpc = 3.0856e22
    condition = np.where((r < r_min) | (r > r_max))[0]
    Theta_r = 3.*((r*Mpc/h0)**2.)/((r_max*Mpc/h0)**3. - (r_min*Mpc/h0)**3.)
    Theta_r[condition] = 0.
    return Theta_r

def get_theta_r_empirical(r, nz, h0):
    """Returns an empirical redshift selection function.

    Parameters
    ----------
    r : array
        Comoving distance.
    nz : float
        redshift selection function.
    h0 : float
        'Little' Hubble constant.
    """
    Mpc = 3.0856e22
    return ((r*Mpc/h0)**2.) * nz / util.integrate(rMpc/h0, ((r*Mpc/h0)**2) * nz, total=True)


def Dz_integrand(a, h0, omega_m, omega_l):
    """Integral for the linear growth function.

    Parameters
    ----------
    a : float
        The scale factor.
    h0 : float
        'Little' Hubble constant.
    omega_m : float
        Present matter density.
    omega_l : float
        Present dark energy density
    """
    return 1./(a*get_H(a, h0, omega_m, omega_l))**3.


def unnorm_Dz(a, h0, omega_m, omega_l):
    """Returns the unnormalised linear growth function.

    Parameters
    ----------
    a : float
        The scale factor.
    h0 : float
        'Little' Hubble constant.
    omega_m : float
        Present matter density.
    omega_l : float
        Present dark energy density
    """
    return get_H(a, h0, omega_m, omega_l)*quad(Dz_integrand, 0., a, args=(h0, omega_m, omega_l))[0]


def get_Dz(z, h0, omega_m, omega_l):
    """Returns the linear growth function.

    Parameters
    ----------
    z : float
        Redshift.
    h0 : float
        'Little' Hubble constant.
    omega_m : float
        Present matter density.
    omega_l : float
        Present dark energy density
    """
    if np.isscalar(z) == True:
        return unnorm_Dz(z2a(z), h0, omega_m, omega_l)/unnorm_Dz(z2a(0.), h0, omega_m, omega_l)
    else:
        Dz = []
        for i in range(0, len(z)):
            Dz.append(unnorm_Dz(z2a(z[i]), h0, omega_m, omega_l))
        Dz = np.array(Dz)
        Dz /= unnorm_Dz(z2a(0.), h0, omega_m, omega_l)
        return Dz


def get_ez(z, omega_m, omega_l):
    """Returns the total amount of matter + dark energy at z.

    Parameters
    ----------
    z : float
        Redshift.
    omega_m : float
        Present matter density.
    omega_l : float
        Present dark energy density
    """
    a = z2a(z)
    return np.sqrt(omega_m*a**-3 + omega_l)


def get_r_integrand(z, omega_m, omega_l):
    """Comoving distance integral.

    Parameters
    ----------
    z : float
        Redshift.
    omega_m : float
        Present matter density.
    omega_l : float
        Present dark energy density
    """
    return 1./get_ez(z, omega_m , omega_l)


def get_r(z, omega_m, omega_l):
    """Returns the comoving distance at z.

    Parameters
    ----------
    z : float
        Redshift.
    omega_m : float
        Present matter density.
    omega_l : float
        Present dark energy density
    """
    if np.isscalar(z) == True:
        return 3000.*quad(get_r_integrand, 0., z, args=(omega_m, omega_l))[0]
    else:
        r = []
        for i in range(0, len(z)):
            r.append(quad(get_r_integrand, 0., z[i], args=(omega_m, omega_l))[0])
        r = np.array(r)
        r *= 3000.
        return r


def get_omega_m_z(z, omega_m, omega_l):
    """Returns the matter density at z.

    Parameters
    ----------
    z : float
        Redshift.
    omega_m : float
        Present matter density.
    omega_l : float
        Present dark energy density
    """
    return (omega_m*(1.+z)**3.)/(omega_m*(1.+z)**3. + omega_l)


def get_fz(z, omega_m, omega_l):
    """Returns the differential of the growth function dlnD/dlna using approximation.

    Parameters
    ----------
    z : float
        Redshift.
    omega_m : float
        Present matter density.
    omega_l : float
        Present dark energy density
    """
    return get_omega_m_z(z, omega_m, omega_l)**0.6
