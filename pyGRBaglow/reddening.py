#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pyGRBaglow import constants as cc
import numpy as np


class reddening(object):
    """
    Class to to model the host and galactic dust and
    gas attenuation of a grb spectrum because
    """

    def __init__(self, wavelength, z, Av):
        """
        Parameters
        ----------
        wavelength: `float` or `array`
        wavelength, in angstroms

        z: `float`
            redshift

        Av : `float`
            Amount of attenuation in V band, in mag
        """

        self.Av = float(Av)
        self.z = float(z)
        self.wavelength = np.atleast_1d(wavelength)

    def NHI_host(self):
        """
        Compute the Hydrogen column density in the host galaxy (in H/cm2)

        Parameters
        ----------
        None

        Returns
        -------
        NHI: `float`
            hydrogen column density, in H/cm2
        """
        # So far we just assume that this is 4 four times higher than
        # in our galaxy.
        # Should be a free parameter if X ray data are used
        return 4. * self.Av * 1.79e21

    def extinction_laws_coeff(self, extinction_law):
        """
        Define the parameters for the extinction law for each cases
        using the "Drude" model by Li and al 2008

        Parameters
        ----------
        extinction_law: `str`
            type of extinction law to use.
            Choices: mw, lmc, smc, linear, calzetti, grb1, grb2

        Returns
        -------
        NHI: `float`
            hydrogen column density, in H/cm2

        """
        if extinction_law.lower() == 'mw':
            c1 = 14.4
            c2 = 6.52
            c3 = 2.04
            c4 = 0.0519
        elif extinction_law.lower() == 'lmc':
            c1 = 4.47
            c2 = 2.39
            c3 = -0.988
            c4 = 0.0221
        elif extinction_law.lower() == 'smc':
            c1 = 38.7
            c2 = 3.83
            c3 = 6.34
            c4 = 0.
        elif extinction_law.lower() == 'linear':
            c1 = 66.2
            c2 = 4.97
            c3 = 22.1
            c4 = 0.
        elif extinction_law.lower() == 'calzetti':
            c1 = 44.9
            c2 = 7.56
            c3 = 61.2
            c4 = 0.
        elif extinction_law.lower() == 'grb1':
            c1 = 0.025
            c2 = 0.048
            c3 = -2.
            c4 = 0.
        elif extinction_law.lower() == 'grb2':
            c1 = 0.015
            c2 = 0.15
            c3 = -2.
            c4 = 0.
        return [c1, c2, c3, c4]

    def Li07(self, extinction_law='smc', Xcut=False):
        """
        Compute the optical depth due to dust reddening

        Parameters
        ----------
        extinction_law : `str`, optional, default: smc
             model of extinction to be used:
             Choices: mw, lmc, smc, linear, calzetti, grb1, grb2

        Xcut: `boolean`, optional, default: False
             Whether to set attenuation to 0 for wavelength below 700 angstrom
             Useful when coupling with X-ray data

        Returns
        ------
        Alambda_over_Av: `float` or `array`
                reddening law

        Trans_dust: `float` or `array`
                Transmission coefficient due to dust reddening
                to be applied to the spectrum
        """
        # Load the coefficient corresponding to the desired
        # extinction curve template
        [c1, c2, c3, c4] = self.extinction_laws_coeff(extinction_law)

        nus = cc.c_light_m_s / (self.wavelength * 1e-10)

        nu = nus * (1 + self.z)
        wavel_mic = self.wavelength * 1e-4 / (1 + self.z)

        # Check if the wavelength is a scalar or an array
        # Otherwise can not iterate over 0-d array
        wavel_mic = np.atleast_1d(wavel_mic)

        # First term accounting for the far UV extinction rise
        UV_extinction = c1 / ((wavel_mic / 0.08)**c2 +
                              (0.08 / wavel_mic)**c2 + c3)

        # Second term accounting for the near-IR/visible extinction
        # Ugly hack to fix non understood indetation pb with pep8
        IR_vis_ext = ((233. * (1. - c1 /
                               (6.88**c2 + 0.145**c2 + c3) - c4 / 4.60))
                      / ((wavel_mic / 0.046)**2. +
                         (0.046 / wavel_mic)**2. + 90.))

        # Third term accounting for the 2175 Angstrom extinction bump
        PAH_bump = c4 / ((wavel_mic / 0.2175)**2. +
                         (0.2175 / wavel_mic)**2. - 1.95)

        # In the rest frame
        Awavel_over_Av = UV_extinction + IR_vis_ext + PAH_bump

        # Set arbitrarily the negative extinction to zero
        w = np.where(Awavel_over_Av < 0.)
        Awavel_over_Av[w] = 0.

        # Applied a cut for wavelength below 700 angstrom
        # Useful when coupling with Xray data
        if Xcut:
            w = np.where(wavel_mic < 0.07)
            Awavel_over_Av[w] = 0.

        # Return optical depth due to dust reddening in funtion of wavelength
        Tau_dust = self.Av / 1.086 * Awavel_over_Av
        Trans_dust = np.exp(-Tau_dust)

        w = np.where(Trans_dust <= 0.)
        Trans_dust[w] = 0.
        w = np.where(Trans_dust > 1.)
        Trans_dust[w] = 1.

        return [Awavel_over_Av, Trans_dust]

    def Pei92(self, Rv='default', ext_law='smc', Xcut=False):
        """
        Extinction laws from Pei 1992 article

        Parameters
        ----------
        Rv: `float`, optional, default: default
            selective attenuation Rv = Av / E(B-V)
            if 'default' set values by default from article
            if a float is given, use this value instead

        ext_law: `str`
            type of extinction law to use.
            Choices: mw, lmc, smc

        Xcut: `boolean`, optional, default: False
             Whether to set attenuation to 0 for wavelength below 700 angstrom
             Useful when coupling with X-ray data

        Returns
        -------
        [Alambda_over_Ab, Trans_dust]

        Alambda_over_Ab : `array`
            atteanuation as a function of wavelength normalise by A_B
            (attenuation in B band)

        Trans_dust: `array`
            transmission through dust as a function of wavelength

        """

        wvl = self.wavelength * 1e-4 / (1 + self.z)
        if ext_law.lower() == 'smc':
            if Rv == 'default':
                Rv = 2.93
            a = [185, 27, 0.005, 0.010, 0.012, 0.03]
            wvl_ = [0.042, 0.08, 0.22, 9.7, 18, 25]
            b = [90, 5.50, -1.95, -1.95, -1.80, 0.0]
            n = [2.0, 4.0, 2.0, 2.0, 2.0, 2.0]

        elif ext_law.lower() == 'lmc':
            if Rv == 'default':
                Rv = 3.16
            a = [175, 19, 0.023, 0.005, 0.006, 0.02]
            wvl_ = [0.046, 0.08, 0.22, 9.7, 18, 25]
            b = [90, 5.5, -1.95, -1.95, -1.8, 0.0]
            n = [2.0, 4.5, 2.0, 2.0, 2.0, 2.0]

        elif ext_law.lower() == 'mw':
            if Rv == 'default':
                Rv = 3.08
            a = [165, 14, 0.045, 0.002, 0.002, 0.012]
            wvl_ = [0.046, 0.08, 0.22, 9.7, 18, 25]
            b = [90, 4.0, -1.95, -1.95, -1.8, 0.0]
            n = [2.0, 6.5, 2.0, 2.0, 2.0, 2.0]

        sums = np.zeros(len(wvl))
        for i in range(len(a)):
            sums += a[i] / ((wvl / wvl_[i])**n[i] +
                            (wvl_[i] / wvl)**n[i] + b[i])

        # Need to check whether extrapolation is needed
        # outside the range defined in Pei92
        Alambda_over_Ab = (1 / Rv + 1) * sums

        # Applied a cut for wavelength below 700 angstrom
        # Useful when coupling with Xray data
        if Xcut:
            w = np.where(wvl < 0.07)
            Alambda_over_Ab[w] = 0

        # Return optical depth due to dust reddening in funtion of wavelength
        Tau_dust = self.Av * Alambda_over_Ab / 1.086

        Trans_dust = np.exp(-Tau_dust)

        Trans_dust[Trans_dust < 0] = 0
        Trans_dust[Trans_dust > 1] = 1

        return [Alambda_over_Ab, Trans_dust]

    def gas_absorption(self, NHx=0.2):
        """
        Compute the optical depth due to gas absorption

        Parameters
        ----------
        NHx: `float`, optional, default: 0.2
             Metal column density from soft Xrays absortpion
             (in units 1e22 /cm2/mag), expressed in units of equivalent
             hydrogen column density assuming solar abundances.
             In Milky Way:
             NHx/Av = 1.7 to 2.2 1e21 cm-2/mag (default set to 2)

        Returns
        ------
        Trans_gas : `array`
             Transmission coefficient due to the gas absorption either
             occuring in our galaxy or within the host.

        """
        nus = cc.c_light_m_s / (self.wavelength * 1e-10)
        Tau_gas = np.zeros(len(nus))

        # Set NHx in units 1e22 /cm2/mag
        NHx *= 1e22

        # Define the equivalent HI column density depending whether we study
        # the host or our galaxy
        """
        if self.z == 0:
             NHI = self.Av * 1.79e21    # cm-2
        else:
             NHI = self.NHI_host()    # cm-2
        """
        NHI = self.Av * NHx

        for i in range(len(nus)):
            # photon frequency (Hz) in the rest frame
            nu = nus[i] * (1 + self.z)
            # photon energy (keV) in the rest frame
            E_kev = nu * cc.H_planck / (1e3 * cc.e_elec)
            E_kev2 = E_kev**2.
            E_kev3 = E_kev**3.

            # if E_kev < 13.6e-3:    #912 A (Lyman limit)
            #     c0=0; c1=0; c2=0
            # 41nm / 410A
            if E_kev < 0.030:
                coeffs = [0, 0, 0, 'H']
            # 12.4 nm
            elif E_kev < 0.100:
                coeffs = [17.3, 608.1, -2150, 'He']
            # 4.37 nm
            elif E_kev < 0.284:
                coeffs = [34.6, 267.9, -476.1, 'C']
            elif E_kev < 0.400:
                coeffs = [78.1, 18.8, 4.3, 'N']
            elif E_kev < 0.532:
                coeffs = [71.4, 66.8, -51.4, 'O']
            elif E_kev < 0.707:
                coeffs = [95.5, 145.8, -61.1, 'Fe-L']
            elif E_kev < 0.867:
                coeffs = [308.9, -380.6, 294.0, 'Ne']
            elif E_kev < 1.303:
                coeffs = [120.6, 169.3, -47.7, 'Mg']
            elif E_kev < 1.840:
                coeffs = [141.3, 146.8, -31.5, 'Si']
            elif E_kev < 2.471:
                coeffs = [202.7, 104.7, -17.0, 'S']
            elif E_kev < 3.210:
                coeffs = [342.7, 18.7, 0.0, 'Ar']
            elif E_kev < 4.038:
                coeffs = [352.2, 18.7, 0.0, 'Ca']
            elif E_kev < 7.111:
                coeffs = [433.9, -2.4, 0.75, 'Fe']
            elif E_kev < 8.331:
                coeffs = [629.0, 30.9, 0.0, 'Ni']
            # 124pm/1.24A
            elif E_kev < 10.:
                coeffs = [701.2, 25.2, 0.0, '...']
            else:
                coeffs = [0., 0., 0., 'None']

            # Figure of M&M
            sige3 = (coeffs[0] + coeffs[1] * E_kev + coeffs[2] * E_kev2)
            # cross section per hydrogen atom /cm2
            sig = sige3 / E_kev3 * 1e-24

            Tau_gas[i] = sig * NHI

        Trans_gas = np.exp(-Tau_gas)

        Trans_gas[Trans_gas < 1e-5] = 0
        Trans_gas[Trans_gas > 1] = 1

        return Trans_gas
