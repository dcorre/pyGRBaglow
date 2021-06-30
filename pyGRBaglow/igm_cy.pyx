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


def meiksin(double[:] wavelength, float redshift, bint lylim=True,
            bint lls_fact=True, bint Xcut=True, float Xlim=10, int num_threads=1):
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

    cdef Py_ssize_t i, j, k
    cdef Py_ssize_t N = len(wavelength)
    cdef Py_ssize_t n_transitions_low = 10
    cdef Py_ssize_t n_transitions_max = 31
    cdef double gamma = 0.2788  # Gamma(0.5,1) i.e., Gamma(2-beta,1) with beta = 1.5
    cdef double n0 = 0.25
    cdef double lambda_limit = 91.2  # Lyman limit in nm
    cdef double term2
    cdef double tau_a, idx, term1
    cdef int[:] w
    cdef double [:] lambda_n = np.zeros(n_transitions_max, dtype=np.float64)
    cdef double [:, ::1] z_n = np.zeros((n_transitions_max, N), dtype=np.float64)
    cdef double[:, ::1] tau_n = np.zeros((n_transitions_max, N), dtype=np.float64)
    cdef double[:] z_l = np.zeros(N, dtype=np.float64)
    cdef double[:] tau_l_igm = np.zeros(N, dtype=np.float64)
    cdef double[:] tau_l_lls = np.zeros(N, dtype=np.float64)
    cdef double[:] term3 = np.zeros(N, dtype=np.float64)
    cdef double[:] term4 = np.zeros(N, dtype=np.float64)
    cdef double[:] tau_taun = np.zeros(N, dtype=np.float64)
    cdef double[:] weight = np.ones(N, dtype=np.float64)
    cdef double[:] igm_transmission = np.zeros(N, dtype=np.float64)

    # From Table 1 in Meiksin (2006), only n >= 3 are relevant.
    # fact has a length equal to n_transitions_low.
    cdef double[:] fact = np.array([1., 1., 1., 0.348, 0.179, 0.109, 0.0722, 0.0508, 0.0373,
                     0.0283])
    if redshift <= 4.:
        tau_a = 0.00211 * pow(1. + redshift,  3.7)
    elif redshift > 4.:
        tau_a = 0.00058 * pow(1. + redshift,  4.5)

    lambda_min_igm = (1.+redshift)*70.

    # LLS Meiksin eq.7
    term1 = gamma - exp(-1.)

    term2 = 0.
    for k in range(n_transitions_low - 1):
        idx = float(k)
        term2 = term2 + pow(-1., idx) / (typed_fact(idx) * (2.*idx - 1.))
       
    for j in prange(N,nogil=True, num_threads=num_threads):
        # Attenuation due to Lyman series lines
        # First, tau_alpha is the mean Lyman alpha transmitted flux,
        # Here n = 2 => tau_2 = tau_alpha
        # Then, tau_n is the mean optical depth value for transitions
        # n = 3 - 9 -> 1
        for i in range(2, n_transitions_max):
            lambda_n[i] = lambda_limit / (1. - 1. / float(i*i))
            z_n[i, j] = (wavelength[j] / lambda_n[i]) - 1.

            if i == 2:
                if redshift <= 4.:
                    tau_n[2, j] = 0.00211 * pow(1. + z_n[2, j], 3.7)
                elif redshift > 4.:
                    tau_n[2, j] = 0.00058 * pow(1. + z_n[2, j], 4.5)
            elif 3 <= i <= 5:
                if z_n[i, j] < 3.:
                    tau_n[i, j] = (tau_a * fact[i] *
                           pow(0.25 * (1. + z_n[i, j]), (1. / 3.)))
                elif z_n[i, j] >= 3:
                    tau_n[i, j] = (tau_a * fact[i] *
                           pow(0.25 * (1. + z_n[i, j]), (1. / 6.)))
            elif 5 < i <= 9:
                tau_n[i,j] = (tau_a * fact[i] *
                            pow(0.25 * (1. + z_n[i, j]), (1. / 3.)))
            elif i > 9:
                tau_n[i, j] = (tau_n[9, j] * 720. /
                            (float(i) * (float(i*i - 1))))
            if z_n[i,j] >= redshift:
                tau_n[i, j] = 0.

            tau_taun[j] = tau_taun[j] + tau_n[i, j]

        # Contribution due to photoelectric absorption: IGM + LLS
        # IGM Meiksin eq.5
        z_l[j] = wavelength[j] / lambda_limit - 1.
        if z_l[j] < redshift:
            tau_l_igm[j] = (0.805 * pow(1. + z_l[j], 3.) *
                            (1. / (1. + z_l[j]) - 1. / (1. + redshift)))
            term3[j] = ((1.+redshift) * pow(wavelength[j]/lambda_limit, 1.5) -
                        pow(wavelength[j] / lambda_limit, 2.5))

        for k in range(n_transitions_low):
            idx = float(k)
            if 1 <= k <= n_transitions_low:
                if z_l[j] < redshift:
                    term4[j] = term4[j] + ((2.* pow(-1., idx) / (typed_fact(idx) * ((6.*idx - 5.)*(2.*idx - 1.)))) *
                                (pow(1.+redshift, 2.5-(3. * idx)) *
                                (wavelength[j] / lambda_limit) ** (3.*idx) -
                                (wavelength[j] / lambda_limit) ** 2.5))

        tau_l_lls[j] = n0 * ((term1 - term2) * term3[j] - term4[j])

        if wavelength[j] < lambda_min_igm:
            if lls_fact == 0:
                # Below lambda_min_igm, transmission set to 0
                weight[j] = 0.
            elif lls_fact == 1:
                weight[j] = pow(wavelength[j]/lambda_min_igm, 2.)
                # Another weight using erf function can be used.
                # However, you would need to add: from scipy.special import erf
                # weight[w] = 0.5*(1.+erf(0.05*(wavelength[w]-lambda_min_igm)))

        igm_transmission[j] = exp(-tau_taun[j]-tau_l_igm[j]-tau_l_lls[j]) * weight[j]

        if igm_transmission[j] < 0.:
            igm_transmission[j]=0.
        if igm_transmission[j] > 1.:
            igm_transmission[j]=1.

        # In case of using this module with data containing XRay data
        # The transmisson has to be set to 1 below an arbitrary wavelength
        # Otherwise there will be no flux in Xray. We chose 10 nm
        if Xcut == 1:
            if wavelength[j] < Xlim * (1.+redshift):
                igm_transmission[j] = 1.

    return igm_transmission
