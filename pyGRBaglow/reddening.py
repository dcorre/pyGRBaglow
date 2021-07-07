#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pyGRBaglow import constants as cc
import numpy as np


def Pei92(wavelength, Av, z, Rv=-99.0, ext_law="smc", Xcut=False):
    """
    Extinction laws from Pei 1992 article

    Parameters
    ----------
    wavelength: `array` or `float`
        wavlength in angstroms

    Av: `float`
        amount of extinction in the V band

    z: `float`
        redshift

    Rv: `float`, optional, default: -99.
        selective attenuation Rv = Av / E(B-V)
        if 'd-99.' set values by default from article
        if a float is given, use this value instead

    ext_law: `str`
        type of extinction law to use.
        Choices: mw, lmc, smc

    Xcut: `boolean`, optional, default: False
         Whether to set attenuation to 0 for wavelength below 700 angstrom
         Useful when coupling with X-ray data

    Returns
    -------
    [Alambda_over_Av, Trans_dust]

    Alambda_over_Av : `array`
        atteanuation as a function of wavelength normalise by Av
        (attenuation in V band)

    Trans_dust: `array`
        transmission through dust as a function of wavelength

    """

    wvl = wavelength * 1e-4 / (1 + z)
    if ext_law.lower() == "smc":
        if Rv == -99.:
            Rv = 2.93
        a = [185, 27, 0.005, 0.010, 0.012, 0.03]
        wvl_ = [0.042, 0.08, 0.22, 9.7, 18, 25]
        b = [90, 5.50, -1.95, -1.95, -1.80, 0.0]
        n = [2.0, 4.0, 2.0, 2.0, 2.0, 2.0]

    elif ext_law.lower() == "lmc":
        if Rv == -99.:
            Rv = 3.16
        a = [175, 19, 0.023, 0.005, 0.006, 0.02]
        wvl_ = [0.046, 0.08, 0.22, 9.7, 18, 25]
        b = [90, 5.5, -1.95, -1.95, -1.8, 0.0]
        n = [2.0, 4.5, 2.0, 2.0, 2.0, 2.0]

    elif ext_law.lower() == "mw":
        if Rv == -99.:
            Rv = 3.08
        a = [165, 14, 0.045, 0.002, 0.002, 0.012]
        wvl_ = [0.046, 0.08, 0.22, 9.7, 18, 25]
        b = [90, 4.0, -1.95, -1.95, -1.8, 0.0]
        n = [2.0, 6.5, 2.0, 2.0, 2.0, 2.0]

    sums = np.zeros(len(wvl))
    for i in range(len(a)):
        sums += a[i] / ((wvl / wvl_[i]) ** n[i] + (wvl_[i] / wvl) ** n[i] + b[i])

    # Need to check whether extrapolation is needed
    # outside the range defined in Pei92
    # convert Alambda_over_Ab to Alambda_over_Av
    Alambda_over_Av = (1.0 / Rv + 1.0) * sums

    # Applied a cut for wavelength below 700 angstrom
    # Useful when coupling with Xray data
    if Xcut:
        w = np.where(wvl < 0.07)
        Alambda_over_Av[w] = 0

    # Return optical depth due to dust reddening in funtion of wavelength
    Tau_dust = Av * Alambda_over_Av / 1.086

    Trans_dust = np.exp(-Tau_dust)

    Trans_dust[Trans_dust < 0] = 0
    Trans_dust[Trans_dust > 1] = 1

    return [Alambda_over_Av, Trans_dust]


def gas_absorption(wavelength, z, NHx=0.2):
    """
    Compute the optical depth due to gas absorption

    Parameters
    ----------
    wavelength: `array` or `float`
        wavlength in angstroms

    z: `float`
        redshift

    NHx: `float`, optional, default: 0.2
         Metal column density from soft Xrays absortpion
         (in units 1e22 cm-2), expressed in units of equivalent
         hydrogen column density assuming solar abundances.
         In Milky Way:
         NHx/Av = 1.7 to 2.2 1e21 cm-2/mag (default set to 2)

    Returns
    ------
    Trans_gas : `array`
         Transmission coefficient due to the gas absorption either
         occuring in our galaxy or within the host.

    """
    nus = cc.c_light_m_s / (wavelength * 1e-10)
    Tau_gas = np.zeros(len(nus))

    for i in range(len(nus)):
        # photon frequency (Hz) in the rest frame
        nu = nus[i] * (1 + z)
        # photon energy (keV) in the rest frame
        E_kev = nu * cc.H_planck / (1e3 * cc.e_elec)
        E_kev2 = E_kev ** 2.0
        E_kev3 = E_kev ** 3.0

        # if E_kev < 13.6e-3:    #912 A (Lyman limit)
        #     c0=0; c1=0; c2=0
        # 41nm / 410A
        if E_kev < 0.030:
            coeffs = [0, 0, 0, "H"]
        # 12.4 nm
        elif E_kev < 0.100:
            coeffs = [17.3, 608.1, -2150, "He"]
        # 4.37 nm
        elif E_kev < 0.284:
            coeffs = [34.6, 267.9, -476.1, "C"]
        elif E_kev < 0.400:
            coeffs = [78.1, 18.8, 4.3, "N"]
        elif E_kev < 0.532:
            coeffs = [71.4, 66.8, -51.4, "O"]
        elif E_kev < 0.707:
            coeffs = [95.5, 145.8, -61.1, "Fe-L"]
        elif E_kev < 0.867:
            coeffs = [308.9, -380.6, 294.0, "Ne"]
        elif E_kev < 1.303:
            coeffs = [120.6, 169.3, -47.7, "Mg"]
        elif E_kev < 1.840:
            coeffs = [141.3, 146.8, -31.5, "Si"]
        elif E_kev < 2.471:
            coeffs = [202.7, 104.7, -17.0, "S"]
        elif E_kev < 3.210:
            coeffs = [342.7, 18.7, 0.0, "Ar"]
        elif E_kev < 4.038:
            coeffs = [352.2, 18.7, 0.0, "Ca"]
        elif E_kev < 7.111:
            coeffs = [433.9, -2.4, 0.75, "Fe"]
        elif E_kev < 8.331:
            coeffs = [629.0, 30.9, 0.0, "Ni"]
        # 124pm/1.24A
        elif E_kev < 10.0:
            coeffs = [701.2, 25.2, 0.0, "..."]
        else:
            coeffs = [0.0, 0.0, 0.0, "None"]

        # Figure of M&M
        sige3 = coeffs[0] + coeffs[1] * E_kev + coeffs[2] * E_kev2
        # cross section per hydrogen atom /cm2
        sig = sige3 / E_kev3  # * 1e-24

        # NHx is given in 1e22 cm-2, and sig should be multiplied by 1e-24
        Tau_gas[i] = sig * NHx * 1e-2

    Trans_gas = np.exp(-Tau_gas)

    Trans_gas[Trans_gas < 1e-5] = 0
    Trans_gas[Trans_gas > 1] = 1

    return Trans_gas
