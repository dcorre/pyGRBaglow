#!/usr/bin/env python
# -*- coding: utf-8 -*-
#cython: boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True

import numpy as np
from pyGRBaglow import constants as cc
from libc.math cimport pow, exp
import cython
from cython.parallel import prange, parallel


cpdef double[:, ::1] Pei92_params(str ext_law="smc"):
    """
    Get parameters for the given extinction law
    Parameters
    ----------
    ext_law: `str`
        type of extinction law to use.
        Choices: mw, lmc, smc
    """
    cdef double[:,::1] params = np.zeros((4,6), dtype=np.float64)

    if ext_law == 'smc':
        params = np.array(
            [[185., 27., 0.005, 0.010, 0.012, 0.03],
             [0.042, 0.08, 0.22, 9.7, 18., 25.],
             [90., 5.50, -1.95, -1.95, -1.80, 0.0],
             [2.0, 4.0, 2.0, 2.0, 2.0, 2.0]],
            dtype=np.float64
        )

    elif ext_law == 'lmc':
        params = np.array(
            [[175., 19., 0.023, 0.005, 0.006, 0.02],
             [0.046, 0.08, 0.22, 9.7, 18., 25.],
             [90., 5.5, -1.95, -1.95, -1.8, 0.0],
             [2.0, 4.5, 2.0, 2.0, 2.0, 2.0]],
            dtype=np.float64
        )

    elif ext_law == 'mw':
        params = np.array(
            [[165., 14., 0.045, 0.002, 0.002, 0.012],
             [0.046, 0.08, 0.22, 9.7, 18., 25.],
             [90., 4.0, -1.95, -1.95, -1.8, 0.0],
             [2.0, 6.5, 2.0, 2.0, 2.0, 2.0]],
            dtype=np.float64
        )

    return params

cpdef double Pei92_float(double wavel, double z, double[:,::1] params,
                         double Rv, str ext_law, bint Xcut) nogil:
    """
    Extinction laws from Pei 1992 article

    Parameters
    ----------
    wavel: `float`
        wavelength in angstroms

    z: `float`
        redshift

    Rv: `float`, optional, default: -99.
        selective attenuation Rv = Av / E(B-V)
        if =-99. set values by default from article
        if a float is given, use this value instead

    ext_law: `str`
        type of extinction law to use.
        Choices: mw, lmc, smc

    Xcut: `boolean`, optional, default: False
         Whether to set attenuation to 0 for wavelength below 700 angstrom
         Useful when coupling with X-ray data

    Returns
    -------
    Alambda_over_Av, Trans_dust

    Alambda_over_Av : `array`
        atteanuation as a function of wavelength normalise by Av
        (attenuation in V band)

    Trans_dust: `array`
        transmission through dust as a function of wavelength

    """
    cdef Py_ssize_t i
    cdef Py_ssize_t N = 4
    cdef double sums = 0.
    cdef double Alambda_over_Av
    cdef double wvl

    wvl = wavel * 1e-4 / (1. + z)
    for i in range(N):
        sums = sums + params[0, i] / (pow(wvl / params[1,i], params[3,i]) +
                                      pow(params[1,i] / wvl, params[3,i]) + params[2,i])
    # Need to check whether extrapolation is needed
    # outside the range defined in Pei92
    # Trnasform Alambda_over_Ab to Alambda_over_Av
    Alambda_over_Av = (1. / Rv + 1.) * sums
    # Assume nothing get absorbed below 700 angstroms as there is the IGM absorption
    # and the analytical formula of Pei92 has a transmission starting increasing at wvl < 800ang
    #if Xcut and wvl < 0.07:
    #    # Applied a cut for wavelength below 700 angstrom
    #    # Useful when coupling with Xray data
    #    Alambda_over_Av = 0.

    return Alambda_over_Av

cpdef double get_trans(double Alambda_over_Av, double Av) nogil:
    cdef double Tau_dust, Trans_dust
    # Return optical depth due to dust reddening in funtion of wavelength
    Tau_dust = Av * Alambda_over_Av / 1.086
    Trans_dust = exp(-Tau_dust)

    if Trans_dust < 0.:
        Trans_dust = 0.
    if Trans_dust > 1.:
        Trans_dust = 1.

    return Trans_dust

def Pei92(double[:] wavelength, double Av, double z, double Rv=-99.,
          str ext_law='smc', bint Xcut=True, int num_threads=1):
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
        if =-99. set values by default from article
        if a float is given, use this value instead

    ext_law: `str`
        type of extinction law to use.
        Choices: mw, lmc, smc

    Xcut: `boolean`, optional, default: False
         Whether to set attenuation to 0 for wavelength below 700 angstrom
         Useful when coupling with X-ray data

    Returns
    -------
    Alambda_over_Av, Trans_dust

    Alambda_over_Av : `array`
        atteanuation as a function of wavelength normalise by Av
        (attenuation in V band)

    Trans_dust: `array`
        transmission through dust as a function of wavelength

    """
    cdef Py_ssize_t i
    cdef Py_ssize_t N = len(wavelength)
    cdef double[:] Alambda_over_Av = np.zeros(N, dtype=np.float64)
    cdef double[:] Trans_dust = np.zeros(N, dtype=np.float64)
    cdef double[:,::1] params #= np.zeros((4,6), dtype=np.float64)

    params = Pei92_params(ext_law)

    if Rv == -99.:
        if ext_law == 'smc':
            Rv = 2.93
        elif ext_law == 'lmc':
            Rv = 3.16
        elif ext_law == 'mw':
            Rv = 3.08

    for i in prange(N, nogil=True, num_threads=num_threads):
        Alambda_over_Av[i] = Pei92_float(
            wavelength[i], z, params, Rv, ext_law, Xcut
        )
        Trans_dust[i] = get_trans(Alambda_over_Av[i], Av)
    
    return Alambda_over_Av, Trans_dust

cpdef double sne_float(double wavel, double z, double[:,::1] params, bint Xcut) nogil:
    """
    Extinction laws from Pei 1992 article

    Parameters
    ----------
    wavel: `float`
        wavelength in angstroms

    z: `float`
        redshift

    Xcut: `boolean`, optional, default: False
         Whether to set attenuation to 0 for wavelength below 700 angstrom
         Useful when coupling with X-ray data

    Returns
    -------
    Alambda_over_Av, Trans_dust

    Alambda_over_Av : `array`
        atteanuation as a function of wavelength normalise by Av
        (attenuation in V band)

    Trans_dust: `array`
        transmission through dust as a function of wavelength

    """
    cdef Py_ssize_t i
    cdef Py_ssize_t N = 4
    cdef double sums = 0.
    cdef double Alambda_over_Av
    cdef double wvl
    cdef double Rv = 2.93
    cdef double scaling_factor = 0.75794464

    wvl = wavel * 1e-4 / (1. + z)
    if wvl <= 1000 * 1e-4:
        for i in range(N):
            sums = sums + params[0, i] / (pow(wvl / params[1,i], params[3,i]) +
                                          pow(params[1,i] / wvl, params[3,i]) + params[2,i])
        # Transform Alambda_over_Ab to Alambda_over_Av
        Alambda_over_Av = (1. / Rv + 1.) * sums * scaling_factor
    else:
        Alambda_over_Av = (
            -2.2113e-05 / pow(wvl, 8)
            + 9.7507e-04 / pow(wvl, 7)
            - 1.7447e-02 / pow(wvl, 6)
            + 1.6186e-01 / pow(wvl, 5)
            - 8.2474e-01 / pow(wvl, 4)
            + 2.262 / pow(wvl, 3)
            - 3.13 / pow(wvl, 2)
            + 2.591 / wvl
            - 6.5916e-01
        )
    # Assume nothing get absorbed below 700 angstroms as there is the IGM absorption
    # and the analytical formula of Pei92 has a transmission starting increasing at wvl < 800ang
    #if Xcut and wvl < 0.07:
    #    # Applied a cut for wavelength below 700 angstrom
    #    # Useful when coupling with Xray data
    #    Alambda_over_Av = 0.

    return Alambda_over_Av

def sne(double[:] wavelength, double Av, double z,
        bint Xcut=True, int num_threads=1):
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

    Xcut: `boolean`, optional, default: False
         Whether to set attenuation to 0 for wavelength below 700 angstrom
         Useful when coupling with X-ray data

    Returns
    -------
    Alambda_over_Av, Trans_dust

    Alambda_over_Av : `array`
        atteanuation as a function of wavelength normalise by Av
        (attenuation in V band)

    Trans_dust: `array`
        transmission through dust as a function of wavelength

    """
    cdef Py_ssize_t i
    cdef Py_ssize_t N = len(wavelength)
    cdef double[:] Alambda_over_Av = np.zeros(N, dtype=np.float64)
    cdef double[:] Trans_dust = np.zeros(N, dtype=np.float64)
    cdef double[:,::1] params #= np.zeros((4,6), dtype=np.float64)

    params = Pei92_params("smc")

    for i in prange(N, nogil=True, num_threads=num_threads):
        Alambda_over_Av[i] = sne_float(wavelength[i], z, params, Xcut)
        Trans_dust[i] = get_trans(Alambda_over_Av[i], Av)

    return Alambda_over_Av, Trans_dust


cpdef (double, double, double) set_gas_absorption_params():
        
    cdef double H_planck = cc.H_planck
    cdef double c_light_m_s = cc.c_light_m_s
    cdef double e_elec = cc.e_elec
    
    return H_planck, c_light_m_s, e_elec

cpdef double gas_absorption_float(double z, double wavel, double NHI, double H_planck,
                                      double c_light_m_s, double e_elec) nogil:
    """
    Compute the optical depth due to gas absorption

    Parameters
    ----------
    wavel: `float`
         wavelength
    
    NHx: `float`, optional, default: 0.2
         Metal column density from soft Xrays absortpion
         (in units 1e22 cm-2), expressed in units of equivalent
         hydrogen column density assuming solar abundances.
         In Milky Way:
         NHx/Av = 1.7 to 2.2 1e21 cm-2/mag (default set to 2)

    Returns
    ------
    Trans_gas : `float`
         Transmission coefficient due to the gas absorption either
         occuring in our galaxy or within the host.

    """

    cdef double freq, coeff_1, coeff_2, coeff_3
    cdef double E_kev, E_kev2, E_kev3, sige3, sig, Tau_gas
    cdef double Trans_gas

    # photon frequency (Hz) in the rest frame
    freq = (1. + z) * c_light_m_s / (wavel * 1e-10)
    # photon energy (keV) in the rest frame
    E_kev = freq * H_planck / (1000. * e_elec)
    E_kev2 = E_kev * E_kev
    E_kev3 = E_kev2 * E_kev

    # Values from Table 2 in http://cdsads.u-strasbg.fr/pdf/1983ApJ...270..119M
    # if E_kev < 13.6e-3:    #912 A (Lyman limit)
    #     c0=0; c1=0; c2=0
    # 41nm / 410A
    if E_kev < 0.030:
        coeff_1 = 0.
        coeff_2 = 0.
        coeff_3 = 0.
        #coeff_4 = 'H'
    # 12.4 nm
    elif E_kev < 0.100:
        coeff_1 = 17.3
        coeff_2 = 608.1
        coeff_3 = -2150.
        #coeff_4 = 'He'
    # 4.37 nm
    elif E_kev < 0.284:
        coeff_1 = 34.6
        coeff_2 = 267.9
        coeff_3 = -476.1
        #coeff_4 = 'C'
    elif E_kev < 0.400:
        coeff_1 = 78.1
        coeff_2 = 18.8
        coeff_3 = 4.3
        #coeff_4 = 'N'
    elif E_kev < 0.532:
        coeff_1 = 71.4
        coeff_2 = 66.8
        coeff_3 = -51.4
        #coeff_4 = 'O'
    elif E_kev < 0.707:
        coeff_1 = 95.5
        coeff_2 = 145.8
        coeff_3 = -61.1
        #coeff_4 = 'Fe-L'
    elif E_kev < 0.867:
        coeff_1 = 308.9
        coeff_2 = -380.6
        coeff_3 = 294.0
        #coeff_4 = 'Ne'
    elif E_kev < 1.303:
        coeff_1 = 120.6
        coeff_2 = 169.3
        coeff_3 = -47.7
        #coeff_4 = 'Mg'
    elif E_kev < 1.840:
        coeff_1 = 141.3
        coeff_2 = 146.8
        coeff_3 = -31.5
        #coeff_4 = 'Si'
    elif E_kev < 2.471:
        coeff_1 = 202.7
        coeff_2 = 104.7
        coeff_3 = -17.0
        #coeff_4 = 'S'
    elif E_kev < 3.210:
        coeff_1 = 342.7
        coeff_2 = 18.7
        coeff_3 = 0.
        #coeff_4 = 'Ar'
    elif E_kev < 4.038:
        coeff_1 = 352.2
        coeff_2 = 18.7
        coeff_3 = 0.
        #coeff_4 = 'Ca'
    elif E_kev < 7.111:
        coeff_1 = 433.9
        coeff_2 = -2.4
        coeff_3 = 0.75
        #coeff_4 = 'Fe'
    elif E_kev < 8.331:
        coeff_1 = 629.0
        coeff_2 = 30.9
        coeff_3 = 0.
        #coeff_4 = 'Ni'
    # 124pm/1.24A
    elif E_kev < 10.:
        coeff_1 = 701.2
        coeff_2 = 25.2
        coeff_3 = 0.0
        #coeff_4 = '...'
    else:
        coeff_1 = 0.
        coeff_2 = 0.
        coeff_3 = 0.
        #coeff_4 = 'None'

    # Figure of M&M
    sige3 = (coeff_1 + coeff_2 * E_kev + coeff_3 * E_kev2)
    # cross section per hydrogen atom 1e-24 cm-2
    sig = sige3 / E_kev3 #* 1e-24

    # NHI given in 1e22 cm-2, and sig should be multiplied by 1e-24
    Tau_gas = sig * NHI * 1e-2

    Trans_gas = exp(-Tau_gas)

    if Trans_gas < 1e-5:
        Trans_gas = 0.
    if Trans_gas > 1.:
        Trans_gas = 1.

    return Trans_gas
    
def gas_absorption(double[:] wavelength, double z, double NHx=0.2, int num_threads=1):
    """
    Compute the optical depth due to gas absorption

    Parameters
    ----------
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
    cdef Py_ssize_t i
    cdef Py_ssize_t N = len(wavelength)
    cdef double[:] Trans_gas = np.zeros(N, dtype=np.float64)
    cdef double H_planck, c_light_m_s, e_elec
    

    H_planck, c_light_m_s, e_elec = set_gas_absorption_params()
            
    for i in prange(N, nogil=True, num_threads=num_threads):
        Trans_gas[i] = gas_absorption_float(z, wavelength[i], NHx, H_planck,
                                                 c_light_m_s, e_elec)

    return Trans_gas
