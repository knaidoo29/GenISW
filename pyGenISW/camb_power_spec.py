import numpy as np
import camb as cb


def get_pk_z_0(h0, omega_b, omega_m, A_s, n_s, npoints=200, minkh=1e-4, maxkh=1e2):
    """Returns the linear power spectra calculated by camb.

    Parameters
    ----------
    h0 : float
        'Little' hubble constant.
    omega_b : float
        Present baryon density.
    omega_m : float
        Present matter density.
    A_s : float
        Amplitude of scalar fluctuations.
    n_s : float
        Spectral tilt.
    npoints : int, optional
        Number of points to calculate the power spectrum.
    minkh : float, optional
        Minimum k to calculate the power spectrum.
    maxkh : float, optional
        Maximum k to calculate the power spectrum.

    Returns
    -------
    kh : array
        Power spectrum k given in units of h/Mpc.
    pk : array
        Power spectrum given in units of (Mpc/h)^3.
    """
    pars = cb.CAMBparams()
    pars.set_cosmology(H0=100.*h0, ombh2=omega_b*h0**2., omch2=(omega_m-omega_b)*h0**2., mnu=0., omk=0.)
    pars.InitPower.set_params(As=A_s, ns=n_s, r=0)
    pars.set_for_lmax(2500, lens_potential_accuracy=0)
    pars.set_matter_power(redshifts=[0.], kmax=10*maxkh)
    #Linear spectra
    pars.NonLinear = cb.model.NonLinear_none
    results = cb.get_results(pars)
    kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=maxkh, npoints=npoints)
    pk = pk.flatten()
    return kh, pk
