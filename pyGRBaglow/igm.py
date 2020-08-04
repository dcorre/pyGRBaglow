#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    This module contains functions used to simulate the InterGalactc
    Medium extinction along the line of sight of a given object such as a GRB.

    meiksin() is taken from the CIGALE code.
"""

from pyGRBaglow import constants as cc
import numpy as np
from scipy.special import factorial


def meiksin(wavelength, redshift, lylim=True,
            lls_fact=False, Xcut=False, Xlim=10):
    """Intergalactic transmission (Meiksin, 2006)

    Compute the intergalactic transmission as described in Meiksin, 2006.

    Parameters
    ----------
    wavelength: `array` or `float`
        wavelength(s) in nm.

    redshift: `float`
        redshift, must be strictly positive.

    lylim: `boolean`, optional, default: True
        Lyman limit on/off. flag whether to include
        photo-electric absorption

    lls_fact: `boolean`, optional, default: False
        factor to increase LLS optical depth: True/False
        If False: set to zero

    Xcut: `boolean`, optional, default: False
          Set transmission to unity below 100nm to allow use of Xray data

    Xlim: `float`, optional, default: 10
          wavelength below which transmission is set to unity (in nm)

    Returns
    -------
    igm_transmission: `array`
        The intergalactic transmission at each input wavelength.

    """
    # Check if the wavelength is a scalar or an array
    # Otherwise can not iterate over 0-d array
    wavelength = np.atleast_1d(wavelength)

    n_transitions_low = 10
    n_transitions_max = 31
    gamma = 0.2788  # Gamma(0.5,1) i.e., Gamma(2-beta,1) with beta = 1.5
    n0 = 0.25
    lambda_limit = 91.2  # Lyman limit in nm

    lambda_n = np.empty(n_transitions_max)
    z_n = np.empty((n_transitions_max, len(wavelength)))
    for n in range(2, n_transitions_max):
        lambda_n[n] = lambda_limit / (1. - 1. / float(n*n))
        z_n[n, :] = (wavelength / lambda_n[n]) - 1.

    # From Table 1 in Meiksin (2006), only n >= 3 are relevant.
    # fact has a length equal to n_transitions_low.
    fact = np.array([1., 1., 1., 0.348, 0.179, 0.109, 0.0722, 0.0508, 0.0373,
                     0.0283])

    # ------------- Attenuation due to Lyman series lines -----------------

    # First, tau_alpha is the mean Lyman alpha transmitted flux,
    # Here n = 2 => tau_2 = tau_alpha
    tau_n = np.zeros((n_transitions_max, len(wavelength)))
    if redshift <= 4:
        tau_a = 0.00211 * np.power(1. + redshift,  3.7)
        tau_n[2, :] = 0.00211 * np.power(1. + z_n[2, :], 3.7)
    elif redshift > 4:
        tau_a = 0.00058 * np.power(1. + redshift,  4.5)
        tau_n[2, :] = 0.00058 * np.power(1. + z_n[2, :], 4.5)

    # Then, tau_n is the mean optical depth value for transitions
    # n = 3 - 9 -> 1
    for n in range(3, n_transitions_max):
        if n <= 5:
            w = np.where(z_n[n, :] < 3)
            tau_n[n, w] = (tau_a * fact[n] *
                           np.power(0.25 * (1. + z_n[n, w]), (1. / 3.)))
            w = np.where(z_n[n, :] >= 3)
            tau_n[n, w] = (tau_a * fact[n] *
                           np.power(0.25 * (1. + z_n[n, w]), (1. / 6.)))
        elif 5 < n <= 9:
            tau_n[n, :] = (tau_a * fact[n] *
                           np.power(0.25 * (1. + z_n[n, :]), (1. / 3.)))
        else:
            tau_n[n, :] = (tau_n[9, :] * 720. /
                           (float(n) * (float(n*n - 1.))))

    for n in range(2, n_transitions_max):
        w = np.where(z_n[n, :] >= redshift)
        tau_n[n, w] = 0.

    # ----------Contribution due to photoelectric absorption: IGM + LLS ------
    tau_l_igm = np.zeros_like(wavelength)
    tau_l_lls = np.zeros_like(wavelength)
    if lylim:
        # IGM Meiksin eq.5
        z_l = wavelength / lambda_limit - 1.
        w = np.where(z_l < redshift)

        tau_l_igm[w] = (0.805 * np.power(1. + z_l[w], 3) *
                        (1. / (1. + z_l[w]) - 1. / (1. + redshift)))

        # LLS Meiksin eq.7
        term1 = gamma - np.exp(-1.)

        n = np.arange(n_transitions_low - 1)
        term2 = np.sum(np.power(-1., n) / (factorial(n) * (2*n - 1)))

        term3 = ((1.+redshift) * np.power(wavelength[w]/lambda_limit, 1.5) -
                 np.power(wavelength[w] / lambda_limit, 2.5))

        term4 = np.sum(np.array(
            [((2.*np.power(-1., n) / (factorial(n) * ((6*n - 5)*(2*n - 1)))) *
              ((1.+redshift) ** (2.5-(3 * n)) *
               (wavelength[w] / lambda_limit) ** (3*n) -
               (wavelength[w] / lambda_limit) ** 2.5))
             for n in np.arange(1, n_transitions_low)]), axis=0)

        tau_l_lls[w] = n0 * ((term1 - term2) * term3 - term4)

        tau_taun = np.sum(tau_n[2:n_transitions_max, :], axis=0)

        lambda_min_igm = (1+redshift)*70.
        w = np.where(wavelength < lambda_min_igm)

    weight = np.ones_like(wavelength)

    if not lls_fact:
        # Below lambda_min_igm, transmission set to 0
        weight[w] = 0
    else:
        weight[w] = np.power(wavelength[w]/lambda_min_igm, 2.)
        # Another weight using erf function can be used.
        # However, you would need to add: from scipy.special import erf
        # weight[w] = 0.5*(1.+erf(0.05*(wavelength[w]-lambda_min_igm)))

    igm_transmission = np.exp(-tau_taun-tau_l_igm-tau_l_lls) * weight

    # In case of using this module with data containing XRay data
    # The transmisson has to be set to 1 below an arbitrary wavelength
    # Otherwise there will be no flux in Xray. We chose 10 nm
    if Xcut:
        w2 = np.where(wavelength < Xlim * (1+redshift))
        igm_transmission[w2] = 1

    # Check for acceptable values
    w = np.where(igm_transmission < 0)
    igm_transmission[w] = 0
    w = np.where(igm_transmission > 1)
    igm_transmission[w] = 1

    return igm_transmission


def madau(wavelength, redshift, lylim=False,
          metals=False, Xcut=False, Xlim=100):
    """
    Compute the intergalactic transmission as described in Madau,
    Madau P. 1995, ApJ, 441, 18.
    For Lyman series line blanketing only Lyalpha (n=2) up to
    Lytheta (n=5) are taken into account.
    Parameters
    ----------
    wavelength : `array`
        wavelength, in obs frame, in angstrom

    redshift : `float`
        redshift of the emitter

    lylim: `boolean`, optional, default: False
        Lyman limit on/off. flag whether to include photo-electric absorption

    metals: `boolean`, optional, default: False
        switch for treating metal blanketing (False = off, True = on)

    Xcut: `boolean`, optional, default: False
        Set transmission to unity below 100nm to allow use of Xray data

    Xlim: `float`, optional, default: 100
        wavelength below which transmission is set to unity (in angstroms)

    Returns
    -------

    Trans_igm : `array`
        Transmission coefficient due to the igm opacity for each wavelength
    """
    # Lyman break
    Lyman_limit = 912.
    a_metal = 0.0017

    Ly_wvl = [1215.67, 1025.72, 972.537, 949.743, 937.803,
              930.748, 926.226, 923.15, 920.963, 919.352,
              918.129, 917.181, 916.429, 915.824, 915.329,
              914.919, 914.576, Lyman_limit]

    Ly_coeff = [0.0036, 0.0017, 0.0011846, 0.0009410, 0.0007960,
                0.0006967, 0.0006236, 0.0005665, 0.0005200, 0.0004817,
                0.0004487, 0.0004200, 0.0003947, 0.000372, 0.000352,
                0.0003334, 0.00031644]

    # initialising tau_eff to zero
    tau_eff = np.zeros(len(wavelength))

    n_transitions_max = len(Ly_wvl)
    # Lyman alpha line blanketing
    tau_n = np.zeros((n_transitions_max, len(wavelength)))

    for n in range(n_transitions_max-1):
        w = np.where((wavelength <= Ly_wvl[n] * (1+redshift))
                     & (wavelength > Ly_wvl[n+1] * (1+redshift)))
        for i in range(n+1):
            tau_n[n, w] += Ly_coeff[i] * (wavelength[w]/Ly_wvl[i])**3.46
            if metals:
                # print ('add metals contribution')
                tau_n[n, w] += a_metal * (wavelength[w]/Ly_wvl[i])**1.68
        tau_eff += tau_n[n, :]

    # ---------Lyman continuum attenuation --------------------
    # approximation given in footnote 3 to the
    # integral in Eq. 16 of Madau (1995).
    # It appears to be a poor approximation for observed wavelengths
    # less than the last transition considered. It should be 912 Angstroms
    # Better to set lylim=False
    w = np.where(wavelength <= Lyman_limit * (1+redshift))
    if lylim:
        xc = wavelength[w] / Lyman_limit
        xm = 1 + redshift
        tau_eff[w] += (0.25 * xc**3. * (xm**0.46 - xc**0.46) + 9.4 * xc**1.5
                       * (xm**0.18 - xc**0.18) - 0.7 * xc**3.
                       * (xc**(-1.32) - xm**(-1.32)) - 0.023
                       * (xm**1.68 - xc**1.68))
    else:
        # transmission set to zero
        tau_eff[w] = np.inf

    # In case of using this module with data containing XRay data
    # The transmisson has to be set to 1 below an arbitrary wavelength
    # Otherwise there will be no flux in Xray. We chose 10 nm
    if Xcut:
        w2 = np.where(wavelength < Xlim * (1 + redshift))
        tau_eff[w2] = 0.

    # Check for acceptable values
    trans_igm = np.exp(-tau_eff)
    w = np.where(trans_igm < 0.)
    trans_igm[w] = 0.
    w = np.where(trans_igm > 1.)
    trans_igm[w] = 1.

    return trans_igm


def dla(wavelength, redshift_dla, NHI):
    """
    Compute the optical depth due to Damped Lyman Alpha systems according
    to Tomani et al 2006 PASJ 58 485-498 June 25.
    Parameters
    ----------
    wavelength: `array`
        wavelength, observe frame, in angstrom

    redshift_dla: `float`
        redshift of the DLA system

    NHI: `float`
        column density of the dla in H/cm2
    Returns
    ------
    Trans_dla: `array`
        Transmission coefficient due to the DLA system
    """
    nus = cc.c_light_m_s / (wavelength * 1e-10)
    f_alpha = 0.4162
    # in cm
    lambda_alpha = 1216e-8
    gu_over_gl = 3.
    # /s
    Lambda_cl_alpha = 1.503e9
    Lambda_alpha = 3. * gu_over_gl**(-1.) * f_alpha * Lambda_cl_alpha
    nu_alpha = cc.c_light_m_s / (1216e-10)

    # photon frequency (Hz) in the rest frame
    nu_rf = nus * (1 + redshift_dla)
    sigma_alpha = ((3. * lambda_alpha**2. * f_alpha * Lambda_cl_alpha
                    / (8.*np.pi)) * Lambda_alpha * (nu_rf / nu_alpha)**4.
                   / (4. * np.pi**2. * (nu_rf - nu_alpha)**2. +
                      Lambda_alpha**2. * (nu_rf / nu_alpha)**(6.) / 4.))
    # ---------------------------------------------------------------------
    # copy right part of sigma_alpha to the left in order to get a symetric
    # cross section. I do not know where the difference comes from.
    # Cross section is not exactly centered on 1216 and depending on the
    # wavelength grids it is more and less symmetric
    w_right = np.where(nu_rf < (nu_alpha))
    sigma_right = sigma_alpha[w_right]
    w_left = np.where(nu_rf > (nu_alpha))
    sigma_left = sigma_alpha[w_left]
    j = min(len(sigma_left), len(sigma_right))
    if len(sigma_left) <= len(sigma_right):
        j = len(sigma_left)
        for i in range(j):
            sigma_left[j-i-1] = sigma_right[i]
        sigma_alpha[w_left] = sigma_left
    else:
        j1 = len(sigma_right)
        j2 = len(sigma_left)
        for i in range(j2):
            if i < j1:
                sigma_left[j2-i-1] = sigma_right[i]
            elif i >= j1:
                # linear interpolation with the last 2 points
                a = ((sigma_left[j2-i+1] - sigma_left[j2-i]) /
                     (nu_rf[j2-i+1] - nu_rf[j2-i]))
                b = sigma_left[j2-i] - a * nu_rf[j2-i]
                sigma_left[j2-i-1] = a * nu_rf[j2-i-1] + b

        sigma_alpha[w_left] = sigma_left
    # ----------------------------------------------------------------

    tau_dla = NHI * sigma_alpha
    # The opacity increases around 200 A due to numerical reasons
    # so below 200 angstrom set to 0
    # w=np.where(nu_rf > cc.c_light_m_s / (200e-10))
    # tau_dla[w] = 0.
    trans_dla = np.exp(-tau_dla)

    w = np.where(trans_dla < 0.)
    trans_dla[w] = 0.
    w = np.where(trans_dla > 1.)
    trans_dla[w] = 1.

    return trans_dla
