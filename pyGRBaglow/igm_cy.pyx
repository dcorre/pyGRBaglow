#!/usr/bin/env python
# -*- coding: utf-8 -*-
#cython: boundscheck=False, wraparound=False, initializedcheck=False, nonecheck=False, cdivision=True

import numpy as np
from libc.math cimport pow, exp
import cython
from cython.parallel import prange, parallel


cdef double typed_fact(double n) nogil:
    """Factorial function for floats"""
    if n <= 1:
        return 1.0
    return n * typed_fact(n-1)


cpdef double[:] set_meiksin_fact():
    # From Table 1 in Meiksin (2006), only n >= 3 are relevant.
    # fact has a length equal to n_transitions_low.
    cdef double[:] fact = np.array([1., 1., 1., 0.348, 0.179, 0.109, 0.0722, 0.0508, 0.0373,
                     0.0283], dtype=np.float64)
    return fact


cpdef double[:] set_meiksin_idx_float(Py_ssize_t n_transitions_max):
    cdef double[:] idx_float=np.arange(n_transitions_max, dtype=np.float64)
    return idx_float


cpdef (double, double, double, double, double, double, double)  set_Meiksin_params(double redshift, double[:] idx_float):
    cdef Py_ssize_t k
    cdef Py_ssize_t n_transitions_low = 10
    cdef double gamma = 0.2788  # Gamma(0.5,1) i.e., Gamma(2-beta,1) with beta = 1.5
    cdef double n0 = 0.25
    cdef double lambda_limit = 91.2  # Lyman limit in nm
    cdef double term2 = 0.
    cdef double tau_a, term1, idx

    if redshift <= 4.:
        tau_a = 0.00211 * pow(1. + redshift,  3.7)
    elif redshift > 4.:
        tau_a = 0.00058 * pow(1. + redshift,  4.5)

    lambda_min_igm = (1.+redshift)*70.

    # LLS Meiksin eq.7
    term1 = gamma - exp(-1.)

    for k in range(n_transitions_low - 1):
        term2 = term2 + pow(-1., idx_float[k]) / (typed_fact(idx_float[k]) * (2.*idx_float[k] - 1.))

    return gamma, n0, lambda_limit, tau_a, lambda_min_igm, term1, term2


cpdef double meiksin_float(double wavelength, double[:] fact, double[:] idx_float,
                          double redshift, double gamma, double n0, double lambda_limit, double tau_a,
                          double lambda_min_igm, double term1, double term2,
                          bint lylim=True, bint lls_fact=False, bint Xcut=True, double Xlim=10.) nogil:
    cdef Py_ssize_t i, k
    cdef Py_ssize_t n_transitions_low = 10
    cdef Py_ssize_t n_transitions_max = 31
    cdef double tau_n = 0.
    cdef double tau_n9 = 0.
    cdef double z_l, lambda_n, z_n
    cdef double tau_l_igm = 0.
    cdef double tau_l_lls =0.
    cdef double term3 = 0.
    cdef double term4 = 0.
    cdef double tau_taun = 0.
    cdef double weight = 1.
    cdef double igm_transmission

    # Attenuation due to Lyman series lines
    # First, tau_alpha is the mean Lyman alpha transmitted flux,
    # Here n = 2 => tau_2 = tau_alpha
    # Then, tau_n is the mean optical depth value for transitions
    # n = 3 - 9 -> 1
    for i in range(2, n_transitions_max):
        lambda_n = lambda_limit / (1. - 1. / (idx_float[i]*idx_float[i]))
        z_n = (wavelength / lambda_n) - 1.

        if i == 2:
            if redshift <= 4.:
                tau_n = 0.00211 * pow(1. + z_n, 3.7)
            elif redshift > 4.:
                tau_n = 0.00058 * pow(1. + z_n, 4.5)
        elif 3 <= i <= 5:
            if z_n < 3.:
                tau_n = (tau_a * fact[i] *
                        pow(0.25 * (1. + z_n), (1. / 3.)))
            elif z_n >= 3:
                tau_n = (tau_a * fact[i] *
                        pow(0.25 * (1. + z_n), (1. / 6.)))
        elif 5 < i <= 9:
            tau_n = (tau_a * fact[i] *
                        pow(0.25 * (1. + z_n), (1. / 3.)))
            if i == 9:
                tau_n9 = tau_n
        elif i > 9:
            tau_n = (tau_n9 * 720. /
                        (idx_float[i] * (idx_float[i] *idx_float[i] - 1.)))
        if z_n >= redshift:
            tau_n = 0.

        tau_taun = tau_taun + tau_n

    # Contribution due to photoelectric absorption: IGM + LLS
    # IGM Meiksin eq.5
    z_l = wavelength / lambda_limit - 1.
    if z_l < redshift:
        tau_l_igm = (0.805 * (1. + z_l) * (1. + z_l) * (1. + z_l) *
                        (1. / (1. + z_l) - 1. / (1. + redshift)))
        term3 = ((1.+redshift) * pow(wavelength/lambda_limit, 1.5) -
                    pow(wavelength / lambda_limit, 2.5))

        for k in range(n_transitions_low):
            if 1 <= k <= n_transitions_low:
                term4 = term4 + ((2.* pow(-1., idx_float[k]) / (typed_fact(idx_float[k]) *
                                                       ((6.*idx_float[k] - 5.)*(2.*idx_float[k] - 1.)))) *
                                 (pow(1.+redshift, 2.5-(3. * idx_float[k])) *
                                  pow(wavelength / lambda_limit, 3.*idx_float[k]) -
                                  pow(wavelength / lambda_limit, 2.5)))

    tau_l_lls = n0 * ((term1 - term2) * term3 - term4)

    if wavelength < lambda_min_igm:
        if lls_fact == 0:
            # Below lambda_min_igm, transmission set to 0
            weight = 0.
        elif lls_fact == 1:
            weight = pow(wavelength/lambda_min_igm, 2.)
            # Another weight using erf function can be used.
            # However, you would need to add: from scipy.special import erf
            # weight[w] = 0.5*(1.+erf(0.05*(wavelength[w]-lambda_min_igm)))

    igm_transmission = exp(-tau_taun-tau_l_igm-tau_l_lls) * weight

    if igm_transmission < 0.:
        igm_transmission=0.
    if igm_transmission > 1.:
        igm_transmission=1.

    # In case of using this module with data containing XRay data
    # The transmisson has to be set to 1 below an arbitrary wavelength
    # Otherwise there will be no flux in Xray. We chose 10 nm
    if Xcut == 1:
        if wavelength < Xlim * (1.+redshift):
            igm_transmission = 1.

    return igm_transmission


def meiksin(double[:] wavelength, double redshift, str unit="nm", bint lylim=True,
            bint lls_fact=False, bint Xcut=True, double Xlim=10., int num_threads=1):
    """Intergalactic transmission (Meiksin, 2006)

    Compute the intergalactic transmission as described in Meiksin, 2006.

    Parameters
    ----------
    wavelength: `array`
        wavelength(s) in nm.

    redshift: `float`
        redshift, must be strictly positive.

    lylim: `boolean`, optional, default: True
        Lyman limit on/off. flag whether to include
        photo-electric absorption. 
        Not used, but kept to remain consistent with arguments of igm.meiksin()

    lls_fact: `boolean`, optional, default: True
        factor to increase LLS optical depth: True/False
        If False: set to zero

    Xcut: `boolean`, optional, default: False
          Set transmission to unity below 100nm to allow use of Xray data

    Xlim: `float`, optional, default: 10
          wavelength below which transmission is set to unity (in nm)

    num_threads: `int`, optional, default: 1
        numbre of threads to use for prange loop

    Returns
    -------
    igm_transmission: `array`
        The intergalactic transmission at each input wavelength.

    """

    cdef Py_ssize_t i
    cdef Py_ssize_t N = len(wavelength)
    cdef Py_ssize_t n_transitions_max = 31
    cdef double[:] igm_transmission = np.zeros(N, dtype=np.float64)
    cdef double gamma, n0, lambda_limit, tau_a, lambda_min_igm, term1, term2
    cdef double corr = 1.
    cdef double[:] fact
    cdef double[:] idx_float 


    fact = set_meiksin_fact()
    idx_float = set_meiksin_idx_float(n_transitions_max)

    gamma, n0, lambda_limit, tau_a, lambda_min_igm, term1, term2 = set_Meiksin_params(redshift, idx_float)

    if unit == "angstroms":
        corr = 1./10.

    for i in prange(N,nogil=True, num_threads=num_threads):
        igm_transmission[i] = meiksin_float(wavelength[i]*corr, fact, idx_float, redshift,
                                           gamma, n0, lambda_limit, tau_a, lambda_min_igm, term1, term2,
                                           lylim=lylim, lls_fact=lls_fact, Xcut=Xcut, Xlim=Xlim)

    return igm_transmission
