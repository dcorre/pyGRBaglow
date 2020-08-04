#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


class Templates(object):
    """ Power law, broken power law..."""

    def __init__(self, F0=500, t0=300, wvl0=6400):
        """
        Parameters
        ----------
        F0 : `float`, optional, default: 500
            flux normaliqation in Jy

        t0: `float`, optional, default: 300
            time corresponding to the normalisation, in s

        wvl0: `float`, optional, default: 6400
            wavelength corresponding to the normalisation, in angstrom

        Returns
        -------
        table : `EventTable`

        """
        self.F0 = F0
        self.t0 = t0
        self.wvl0 = wvl0

        return None

    def SPL(self, wvl, t, alpha, beta):
        """ Simple Power Law.

        Parameters
        ----------
        wvl : `array` or `float`
            wavelength at which to compute the flux
            Must be same unit as wwl0

        t: `array` or `float`
            time at which to compute the flux
            Same unit as t0

        alpha: `float`
            temporal index

        beta: `float`
            spectral index

        Returns
        -------
        flux: `array` or `float`

        """

        flux = self.F0 * (t / self.t0)**(-alpha) * (wvl / self.wvl0)**(beta)
        return flux

    def BPL(self, wvl, t, alpha1, alpha2, beta, s):
        """ Broken Power Law.

        Parameters
        ----------
        wvl : `array` or `float`
            wavelength at which to compute the flux
            Must be same unit as wwl0

        t: `array` or `float`
            time at which to compute the flux
            Same unit as t0

        alpha1: `float`
            temporal index of first power law

        alpha2: `float`
            temporal index of second power law

        beta: `float`
            spectral index

        s: `float`
            smoothing parameter connecting the two powerlaws

        Returns
        -------
        flux: `array` or `float`

        """
        F = self.F0 * (wvl / self.wvl0)**(beta) \
            * ((t / self.t0)**(-s * alpha1)
               + (t / self.t0)**(-s * alpha2))**(-1/s)
        return F

    def light_curve(self, wavelength, time, params, model='SPL'):
        """ Compute light curves

        Parameters
        ----------
        wavelength : `array` or `float`
            wavelength at which to compute the flux
            Must be same unit as wwl0

        time: `array` or `float`
            time at which to compute the flux
            Same unit as t0

        params: `list`
            contains list of params for te given model (SPL/BPL)

        model: `str`, optional, default: SPL
            model to select between SPL and BPL

        Returns
        -------
        flux: `array` or `float`
            if wvl and t are float it returns a float
            otherwise it is a 2D array containing the flux at each given
            wavelength and time

        """
        time_series = np.atleast_1d(time)
        wavelength = np.atleast_1d(wavelength)
        lc = []
        t1 = len(time)
        wvl1 = len(wavelength)
        for t in range(t1):
            sed = []
            for wvl in range(wvl1):
                if model == 'SPL':
                    alpha, beta = params
                    SED = self.SPL(wavelength[wvl], time[t], alpha, beta)

                elif model == 'BPL':
                    alpha1, alpha2, beta, s = params
                    SED = self.BPL(wavelength[wvl], time[t],
                                   alpha1, alpha2, beta, s)
                # Flux in obs frame in same unit as F0
                sed.append(SED)
            # Flux in obs frame in same unit as F0
            lc.append(sed)
        return np.array(lc)
