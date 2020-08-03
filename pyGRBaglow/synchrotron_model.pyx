#!/usr/bin/env python
# -*- coding: utf-8 -*-
#cython: boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True

import numpy as np
from astropy.cosmology import Planck15 as cosmo
from pyGRBaglow import constants as cc
from libc.math cimport pow, exp, pi
from libc.stdio cimport printf
import cython
from cython.parallel import prange, parallel

cdef class fireball_afterglow(object):
    """  Class to to model a grb afterglow according to the fireball model """
    cdef long double E_iso, Msun, Mpc_to_cm, A
    cdef double Mdot_loss, DL
    cdef float n0, eps_b, eps_e, eta, p, Y, Vw, z, td, nu
    cdef int ism_type, disp, num_threads

    def __cinit__(self, float n0=0.001, float eps_b=0.001, float eps_e=0.001,
                 long double E_iso=2*pow(10., 52.), float eta=0.1, float p=2.4,
                 float Y=0., double Mdot_loss=pow(10,-5), float Vw=1000.,
                 float z=1., float td=1., float nu=pow(10.,15.),
                 int ism_type=0, int disp=0, int num_threads=1):
        self.n0 = n0
        self.eps_b = eps_b
        self.eps_e = eps_e
        self.E_iso = E_iso
        self.eta = eta
        self.p = p
        self.Y = Y
        self.Mdot_loss = Mdot_loss
        self.Vw = Vw
        self.z = z
        self.td = td
        self.nu = nu
        self.ism_type = ism_type
        self.disp = disp
        self.num_threads = num_threads
        
        self.Msun = cc.m_sun
        self.Mpc_to_cm = cc.Mpc_cm
        
        self.A = (self.Mdot_loss * self.Msun * pow(10.,3.)/ 365.25 / 86400.) \
             / (4. * pi * self.Vw * pow(10.,5.))

    cdef float granot_sari(self, float td=1., float nu=1*pow(10.,15.)) nogil:
        """ From Granot et Sari """

        cdef long double E52, d, A_star, nu14, eps_ebar, evolution_param
        cdef long double nu1, nu2, nu3, nu4, nu5, nu6, nu7, nu8, nu9, nu10, nu11
        cdef long double f1, f2, f3, f4, f5, f6, f7 , f8, f9, f10, f11
        cdef long double Ftilde1, Ftilde2, Ftilde3, Ftilde4, Ftilde5, Ftilde6
        cdef long double Ftilde7 , Ftilde8, Ftilde9, Ftilde10, Ftilde11
        cdef long double F1, F2, F3, F4, F5, F6, F7
        cdef long double phi1, phi2, phi3, phi4
        cdef long double nuac, nusa, nuc, num, f, ff, fsa, fm, fc, fac
        cdef long double nub, fext, bet1, bet2, s
        cdef double t_obs, 
        
        cdef float s1, s2, s3, s4, s5, s6, s7 , s8, s9, s10, s11
        cdef float bet1_1, bet1_2, bet1_3, bet1_4, bet1_5, bet1_6, bet1_7
        cdef float bet1_8, bet1_9, bet1_10, bet1_11
        cdef float bet2_1, bet2_2, bet2_3, bet2_4, bet2_5, bet2_6, bet2_7
        cdef float bet2_8, bet2_9, bet2_10, bet2_11
         
        cdef int spectrum_case, evolution_case

        #cdef long double Mpc_to_cm = 3.085677e+24
        
        t_obs = td * 86400.  # duration in seconds

        E52 = (1-self.eta)/self.eta * self.E_iso*1*pow(10.,-52.)
        #d = DL*cc.Mpc_cm / (1e28)
        # Divide by 2e28 to be consistent with Turpin's catalogue.
        d = self.DL * self.Mpc_to_cm / pow(10.,28.)
        #d = 4.9190e+28 / (2e28)

        eps_ebar = self.eps_e *(self.p-2.)/(self.p-1.)

        nu14 = nu / pow(10., 14.)

        # Constant ISM
        if self.ism_type == 0:
   
            if self.disp == 1:
                 printf("CONSTANT ISM: \n")
                 printf('Forward shock, at tobs = %f seconds\n',t_obs)
   
            # Granot et Sari 2002 (Section 5, page 825)
            evolution_param = self.n0 * E52**(4./7) * self.eps_b**(9./7)
            if evolution_param < 18:
                 evolution_case = 512
            else:
                 evolution_case = 432
   
            # Granot et Sari 2002 (Table 2)
    
            nu1 = 1.24 * (pow(self.p-1,3./5.) / pow(3.*self.p+2,3./5.)) * pow(10.,9.) \
                  * pow(1.+self.z,-1.) * pow(eps_ebar,-1.) * pow(self.eps_b, 1./5.) \
                  * pow(self.n0, 3./5.) * pow(E52, 1./5.)

            nu2 = 3.73 * (self.p-0.67) * pow(10., 15.) * pow(1.+self.z, 0.5) \
                  * pow(E52, 0.5) * pow(eps_ebar, 2.) * pow(self.eps_b, 0.5) \
                  * pow(td, -3./2.)
            
            nu3 = 6.37 * (self.p-0.46) * pow(10., 13.) * exp(-1.16 * self.p) \
                  * pow(1.+self.z, -0.5) * pow(self.eps_b, -3./2.) * pow(self.n0, -1.) \
                  * pow(E52, -0.5) * pow(td, -0.5)
            
            nu4 = 5.04 * (self.p-1.22) * pow(10., 16.) * pow(1.+self.z, 0.5) \
                  * pow(eps_ebar, 2.) * pow(self.eps_b, 0.5) * pow(E52, 0.5) \
                  * pow(td,-3./2.)

            nu5 = 3.59 * (4.03-self.p) * pow(10., 9.) * exp(2.34*self.p) \
                  * pow( pow(eps_ebar, 4.*(self.p-1.)) * pow(self.eps_b, self.p+2.) * 
                       pow(self.n0, 4.) * pow(E52, self.p+2) / pow(1.+self.z, (6.-self.p))
                       / pow(td, 3.*self.p+2.) , 1./(2.*(self.p+4.)))
            
            nu6 = 3.23 * (self.p-1.76) * pow(10.,12.) \
                  * pow( pow(eps_ebar, 4.*(self.p-1.)) * pow(self.eps_b, self.p-1.) \
                         * pow(self.n0, 2.) * pow(E52, self.p+1) / pow(1.+self.z, 7.-self.p) 
                         / pow(td, 3.*(self.p+1.)), 1./(2.*(self.p+5.)))
            
            nu7 = 1.12 * pow(3.*self.p-1., 8./5.) / pow(3.*self.p+2,8./5) \
                  * pow(10.,8.) * pow(1.+self.z, -13./10) * pow(eps_ebar,-8./5.) \
                  * pow(self.eps_b, -2./5.) * pow(self.n0,3./10) * pow(E52, -1./10.) \
                  * pow(td, 3./10)
            
            nu8 = 1.98 * pow(10.,11.) * pow(1.+self.z,-0.5) * pow(self.n0, 1./6.) \
                  * pow(E52, 1./6.) * pow(td, -0.5)
            
            nu9 = 3.94 * (self.p-0.74) * pow(10., 15.) * pow(1.+self.z, 0.5) \
                  * pow(eps_ebar,2.) * pow(self.eps_b,0.5) * pow(E52,0.5) \
                  * pow(td, -3./2.)
            
            # Following line must be updated
            nu10 = 1.32 * pow(10.,10.) * pow(1.+self.z, -0.5) * pow(self.eps_b, 6./5.) \
                   * pow(self.n0, 11./10.) * pow(E52, 7./10.) * pow(td, -0.5)
            # Following line must be updated
            nu11 = 5.86 * pow(10., 12.) * pow(1.+self.z, -0.5) * pow(self.eps_b, -3./2.) \
                   * pow(self.n0,-1.) * pow(E52, -0.5) * pow(td, -0.5)
            
            f1 = 0.647 * pow(self.p-1., 6./5.) / ( (3.*self.p-1.) * pow(3.*self.p+2,1./5.) )\
                 * pow(1.+self.z, 0.5) * pow(eps_ebar,-1.) * pow(self.eps_b, 2./5.) \
                 * pow(self.n0, 7./10.) * pow(E52, 9./10.) * pow(td, 0.5) * pow(d, -2.)
            
            f2 = 9.93 * (self.p+0.14) * (1.+self.z) * pow(self.eps_b, 0.5) \
                * pow(self.n0,0.5) * E52 * pow(d,-2.)
                       
            f3 = 4.68 * exp(4.82 * (self.p-2.5)) * pow(10.,3.) * pow(1.+self.z, (self.p+1.)/2.) \
                 * pow(eps_ebar, self.p-1) * pow(self.eps_b, self.p-0.5) \
                 * pow(self.n0, self.p/2.) * pow(E52, (self.p+1.)/2.) * pow(td,(1.-self.p)/2.) \
                 * pow(d,-2.)
            
            f4 = 3.72 * (self.p-1.79) * pow(10., 15.) * pow(1.+self.z,7./2.) \
                 * pow(eps_ebar, 5.) * self.eps_b * pow(self.n0,-0.5) \
                 * pow(E52,3./2.) * pow(td,-5./2) * pow(d,-2.)
                       
            f5 = 20.8 * (self.p-1.53) * exp(2.56*self.p) * pow(d,-2.) \
                 * pow( pow(1.+self.z, 7.*self.p+3.) * pow(self.eps_b, 2.*self.p+3.) \
                        * pow(E52, 3.*self.p+7.) / pow(eps_ebar, 10.*(1.-self.p)) \
                        / pow(td, 5.*(self.p-1.)), 1./(2.*(self.p+4.)))
            
            f6 = 76.9 * (self.p-1.08) * exp(2.06*self.p) *pow(d,-2.) \
                 * pow( pow(1.+self.z, 7.*self.p+5.) * pow(self.eps_b, 2.*self.p-5.) \
                        * pow(E52, 3.*self.p+5.) / pow(eps_ebar, 10.*(1.-self.p)) \
                        / pow(self.n0, self.p) / pow(td, 5.*(self.p-1)), 1./(2.*(self.p+5)))
            
            f7 = 5.27 * pow(3.*self.p-1., 11./5.) / pow(3.*self.p+2., 11./5.) \
                 * pow(10., -3.) * pow(1.+self.z,-1./10.) * pow(eps_ebar,-11./5.) \
                 * pow(self.eps_b, -4./5.) * pow(self.n0, 1./10.) * pow(E52,3./10.) \
                 * pow(td,11./10.) * pow(d,-2.)
            
            f8 = 154. * (1.+self.z) * pow(self.eps_b,-1./4) * pow(self.n0, -1./12.) \
                 * pow(E52,2./3.) * pow(d,-2.)
                       
            f9 = 0.221 * (6.27-self.p) * pow(1.+self.z,0.5) * pow(eps_ebar,-1.) \
                 * pow(self.eps_b,-0.5) * pow(E52,0.5) * pow(td, 0.5) * pow(d,-2.)
                       
            f10 = 3.72 * (1.+self.z) * pow(self.eps_b, 7./5.) * pow(self.n0, 6./5.) \
                  * pow(E52, 7./5) * pow(d,-2.)
                       
            f11 = 28.4 * (1.+self.z) * pow(self.eps_b, 0.5) * pow(self.n0, 0.5) \
                  * E52 * pow(d,-2.)  
   
            #print (hex(id(eps_ebar)))
            s1 = 1.64
            s2 = 1.84 - 0.40 * self.p
            s3 = 1.15 - 0.06 * self.p
            s4 = 3.44 * self.p - 1.41
            s5 = 1.47 - 0.21 * self.p
            s6 = 0.94 - 0.14 * self.p
            s7 = 1.99 - 0.04 * self.p
            s8 = 0.907
            s9 = 3.34 - 0.82 * self.p
            s10 = 1.213
            s11 = 0.597
   
            # Compute the choise of the spectrum case
            spectrum_case = 0
            if evolution_case == 512:
                if nu11 <= nu9:
                    spectrum_case = 5
                else:
                    if nu1 <= nu2:
                        spectrum_case = 1
                    else:
                        spectrum_case = 2
   
            elif evolution_case == 432:
                if nu8 <= nu9:
                     spectrum_case = 4
                else:
                    if nu5 <= nu3:
                        spectrum_case = 2
                    else:
                        spectrum_case = 3
        else:
            # Massive star progenitor surounded by its preexplosion Wind
            if self.disp:
                printf("WINDY environment: \n")
                printf('Forward shock, at tobs = %f seconds\n',t_obs)
   
            A_star = self.A /(5. * pow(10.,11.))
   
            evolution_param = A_star * pow(eps_ebar,-1.) * pow(E52,-3./7) * pow(self.eps_b,2./7)
   
            if evolution_param > 100.:
                 evolution_case = 4512
            else:
                 evolution_case = 432

            # Granot et Sari 2002 (Table 2)
            nu1 = 8.31 * pow(self.p-1., 3./5.) / pow(3.*self.p+2,3./5.) * pow(10.,9.) \
                  * pow(1.+self.z,-2./5.) * pow(eps_ebar,-1.) * pow(self.eps_b,1./5.) \
                  * pow(A_star,6./5.) * pow(E52,-2./5.) * pow(td,-3./5.)
           
            nu2 = 4.02 * (self.p-0.69) * pow(10.,15.) * pow(1.+self.z, 0.5) \
                  * pow(E52, 0.5) * pow(eps_ebar, 2.) * pow(self.eps_b, 0.5) * pow(td, -3./2.)
            
            nu3 = 4.40 * (3.45-self.p) * pow(10.,10.) * exp(0.45*self.p) * pow(1.+self.z, -3./2.) \
                  * pow(self.eps_b, -3./2.) * pow(A_star,-2.) * pow(E52, 0.5) * pow(td,0.5)
            
            nu4 = 8.08 * (self.p-1.22) * pow(10.,16.) * pow(1.+self.z, 0.5) \
                  * pow(eps_ebar,2.) * pow(self.eps_b, 0.5) * pow(E52, 0.5) * pow(td,-3./2.)
                       
            nu5 = 1.58 * (4.10-self.p) * pow(10., 10.) * exp(2.16*self.p) \
                  * pow( pow(eps_ebar,4.*(self.p-1.)) * pow(self.eps_b, self.p+2.) \
                         * pow(A_star, 8.) / pow(1.+self.z, 2.-self.p) / pow(E52,2.-self.p) \
                         / pow(td, 3.*(self.p+2.)), 1./(2.*(self.p+4.)))

            nu6 = 4.51 * (self.p-1.73) * pow(10.,12.) \
                  * pow( pow(eps_ebar, 4.*(self.p-1.)) * pow(self.eps_b, self.p-1.) \
                         * pow(A_star, 4.) * pow(E52, self.p-1.) / pow(1.+self.z,5.-self.p) \
                         / pow(td, 3.*self.p+5.), 1./(2.*(self.p+5.)))

            nu7 = 1.68 * pow(3.*self.p-1., 8./5.) / pow(3.*self.p+2., 8./5) * pow(10.,8.) \
                  * pow(1.+self.z,-1.) * pow(eps_ebar,-8./5) * pow(self.eps_b,-2./5.) \
                  * pow(A_star,3./5.) * pow(E52,-2./5.)

            nu8 = 3.15 * pow(10.,11.) * pow(1.+self.z,-1./3.) * pow(A_star,1./3.) \
                  * pow(td,-2./3.)
            
            nu9 = 3.52 * (self.p-0.31) * pow(10.,15.) * pow(1.+self.z,0.5) \
                  * pow(eps_ebar,2.) * pow(self.eps_b,0.5) * pow(E52, 0.5) * pow(td,-3./2)
            
            nu10 = 1.32 * pow(10.,10.) * pow(1.+self.z, -0.5) * pow(self.eps_b,6./5.) \
                   * pow(self.n0, 11./10.) * pow(E52, 7./10.) * pow(td,-0.5)
            
            nu11 = 5.86 * pow(10.,12.) * pow(1.+self.z, -0.5) * pow(self.eps_b,-3./2.) \
                   * pow(self.n0,-1.) * pow(E52, -0.5) * pow(td, -0.5)

            f1 = 9.19 * pow(self.p-1,6./5.) / ( (3.*self.p-1.) * pow(3*self.p+2,1./5.) ) \
                 * pow(1.+self.z, 6./5.) * pow(eps_ebar,-1.) * pow(self.eps_b, 2./.5) \
                 * pow(A_star, 7./5.) * pow(E52,1./5.) * pow(td,-1./5.) * pow(d,-2.)
            
            f2 = 76.9 * (self.p + 0.12) * pow(1.+self.z,3./2.) * pow(self.eps_b,1./2.) \
                 * A_star * pow(E52, 0.5) * pow(td, -0.5) * pow(d, -2.)
   
            f3 = 8.02 * exp(7.02 * (self.p-2.5)) *pow(10.,5.) * pow(1.+self.z, self.p+1./2.) \
                 * pow(eps_ebar,self.p-1) * pow(self.eps_b,self.p-1./2.) * pow(A_star,self.p) \
                 * pow(E52, 0.5) * pow(td,0.5-self.p) * pow(d,-2.)
            
            f4 = 3.04 * (self.p-1.79) * pow(10.,15.) * pow(1.+self.z,3.) * pow(eps_ebar,5.) \
                 * self.eps_b * pow(A_star,-1.) * pow(E52, 2.) * pow(td,-2.) * pow(d,-2.)
            
            f5 = 158. * (self.p-1.48) * exp(2.24 * self.p) * pow(d,-2.) \
                 * pow( pow(1.+self.z, 6.*self.p+9) * pow(self.eps_b,2.*self.p+3.) \
                        * pow(E52,4.*self.p+1.) / pow(eps_ebar, 10.*(1.-self.p)) \
                        / pow(A_star, 2.*(self.p-6.)) * pow(td, 4*self.p+1.), 1./(2.*(self.p+4.)))

            f6 = 8.6 * (self.p-1.12) * exp(1.89 * self.p) * pow(d,-2.) \
                 * pow( pow(1.+self.z,6.*self.p+5.) *pow(self.eps_b,2.*self.p-5.) \
                        * pow(E52, 4.*self.p+5.) / pow(eps_ebar, 10.*(1.-self.p)) \
                        / pow(A_star, 2.*self.p) / pow(td,4.*self.p-5.), 1./(2.*(self.p+5.)))

            f7 = 3.76 * pow(3.*self.p-1,11./5.) / (pow(3.*self.p+2.,11./5.) ) * pow(10.,-3.) \
                 * pow(eps_ebar,-11./5.) * pow(self.eps_b, -4./5.) * pow(A_star,1./5.) \
                 * pow(E52,1./5.) * td * pow(d, -2.)

            f8 = 119. * pow(1.+self.z, 11./12.) * pow(self.eps_b,-1./4.) * pow(A_star,-1./6.) \
                 * pow(E52,3./4.) * pow(td,1./12.) * pow(d,-2.)

            f9 = 0.165 * (7.14-self.p) * pow(1.+self.z,0.5) * pow(eps_ebar,-1.) \
                 * pow(self.eps_b,-0.5) * pow(E52,0.5) * pow(td,0.5) * pow(d,-2.)
            
            # Following line must be updated# Following line must be updated
            f10 = 3.72 * (1.+self.z) * pow(self.eps_b,7./5.) * pow(self.n0,6./5) \
                  * pow(E52, 7./5.) * pow(d,-2.)

            # Following line must be updated
            f11 = 28.4 * (1.+self.z) * pow(self.eps_b,0.5)* pow(self.n0,0.5) * E52 * pow(d,-2.)
   
            s1 = 1.06
            s2 = 1.76 - 0.38 * self.p
            s3 = 0.80 - 0.03 * self.p
            s4 = 3.63 * self.p - 1.60
            s5 = 1.25 - 0.18 * self.p
            s6 = 1.04 - 0.16 * self.p
            s7 = 1.97 - 0.04 * self.p
            s8 = 0.893
            s9 = 3.68 - 0.89 * self.p
            # Following line must be updated
            s10 = 1.213
            # Following line must be updated
            s11 = 0.597
   
            #Compute the choise of the spectrum case
            spectrum_case = 0
   
            if evolution_case == 4512:
                if nu11 <= nu8:
                    spectrum_case = 4
                else:
                    if nu11 <= nu9:
                        spectrum_case = 5
                    else:
                        if nu1 <= nu2:
                            spectrum_case = 1
                        else:
                            spectrum_case = 2
   
            elif evolution_case == 432:
                if nu8 <= nu9:
                    spectrum_case = 4
                else:
                    if nu5 <= nu3:
                        spectrum_case = 2
                    else:
                        spectrum_case = 3

        # Common part for ISM and Wind  
        if self.disp:
            printf('Evolution case : %d \n', evolution_case)
            printf('Spectrum Type : %d \n', spectrum_case)
   
        # Inverse Compton discussed in page 827 
        # Y = eps_e /eps_b if eps_e << eps_b
        # Y = sqrt(eps_e/eps_b) if eps_e >> eps_b (Sari, Piran et Narayan 1996)
   
        if self.Y > 0. and self.eps_b < 10.* self.eps_e:
            if self.disp:
                printf('Inverse Compton occurs')
   
            nu1 = nu1 * (1. + self.Y)
            nu2 = nu2 * (1. + self.Y)
            nu3 = nu3 * (1. + self.Y)
            nu4 = nu4 * (1. + self.Y)
            nu5 = nu5 * (1. + self.Y)
            nu6 = nu6 * (1. + self.Y)
            nu7 = nu7 * (1. + self.Y)
            nu8 = nu8 * (1. + self.Y)
            nu9 = nu9 * (1. + self.Y)
            nu10 = nu10 * (1. + self.Y)
            nu11 = nu11 * (1. + self.Y)
   
        else:
            if self.disp:
                printf('No inverse Compton expected')
   
        # Granot et SAri 2002 (Table 2)
        bet1_1 = 2.
        bet1_2 = 1. / 3.
        bet1_3 = (1 - self.p) / 2.
        bet1_4 = 2.
        bet1_5 = 5. / 2.
        bet1_6 = 5. / 2.
        bet1_7 = 2.
        bet1_8 = 11. / 8.
        bet1_9 = -1. / 2.
        bet1_10 = 11. / 8.
        bet1_11 = 1. / 3.
    
        bet2_1 = 1. / 3.
        bet2_2 = (1. - self.p) / 2.
        bet2_3 = -self.p / 2.
        bet2_4 = 5. / 2.
        bet2_5 = (1. - self.p) / 2.
        bet2_6 = -self.p / 2.
        bet2_7 = 11. / 8.
        bet2_8 = -1. / 2.
        bet2_9 = -self.p / 2.
        bet2_10 = 1. / 3.
        bet2_11 = -1. / 2.
   
        nuac = 0.
        nusa = 0.
        nuc = 0.
        num = 0.
        f = 0.
        ff = 0.

        # Validity of Spectrum 1
        if spectrum_case == 1:
            nusa = nu1
            fsa = f1
            num = nu2
            fm = f2
            nuc = nu3
            fc = f3
            # Prescription for the broad band spectra (formulae 4-9)
            Ftilde2 = pow(1. + pow(nu / nu2, s2 * (bet1_2 - bet2_2)), -1. / s2 )
            Ftilde3 = pow(1. + pow(nu / nu3, s3 * (bet1_3 - bet2_3)), -1. / s3 )
            nub = nu1
            fext = f1
            bet1 = bet1_1
            bet2 = bet2_1
            s = s1
            F1 = fext *pow( pow(nu/nub, -s*bet1) + pow(nu/nub, -s*bet2), -1./s)
            ff = F1 * Ftilde2 * Ftilde3

        # Validity of spectrum 2
        if spectrum_case == 2:
            num = nu4
            fm = f4
            nusa = nu5
            fsa = f5
            nuc = nu3
            fc = f3
            # Prescription for the broadband spectra (formulae 4-9)
            Ftilde5 = pow(1. + pow(nu/nu5, s5 *(bet1_5 - bet2_5)), -1./s5 )
            Ftilde3 = pow(1. + pow(nu/nu3, s3 *(bet1_3 - bet2_3)), -1./s3 )
            nub = nu4
            fext = f4
            bet1 = bet1_4
            bet2 = bet2_4
            s = s4
            if nub == nu4:
                phi4 = (nu/nu4) # (3)
                F4 = f4 * (phi4 * phi4 * exp(-s * pow(phi4,2./3.)) + pow(phi4, 5./2.))
            else:
                F4 = fext * pow( pow(nu/nub, -s*bet1) + pow(nu/nub, -s*bet2), -1./s)
            ff = F4 * Ftilde5 * Ftilde3

        # Validity of spectrum 3
        if spectrum_case == 3:
            num = nu4
            fm = f4
            nusa = nu6
            fsa = f6
            # Prescription for the broadband spectra (formulae 4-9)
            Ftilde6 = pow(1. + pow(nu/nu6, s6 *(bet1_6 -bet2_6 )), -1./s6 )
            nub = nu4
            fext = f4
            bet1 = bet1_4
            bet2 = bet2_4
            s = s4
            if nub == nu4:
                phi4 = (nu/nu4) # (3)
                F4 = f4 * (phi4  * phi4 * exp(-s * pow(phi4,2./3.)) + pow(phi4,5./2.))
            else:
                F4 = fext * pow( pow(nu/nub,-s*bet1) + pow(nu/nub,-s*bet2), -1./s)
            ff = F4 * Ftilde6

        # Validity of spetrum 4 
        if spectrum_case == 4:
            nuac = nu7
            fac = f7
            nusa = nu8
            fsa = f8
            num = nu9
            fm = f9
            # Prescription for the broadband spectra (formulae 4-9)
            Ftilde8 = pow(1. + pow(nu/nu8, s8 *(bet1_8 -bet2_8)), -1./s8 )
            Ftilde9 = pow(1. + pow(nu/nu9, s9 *(bet1_9 -bet2_9)), -1./s9 )
            nub = nu7
            fext = f7
            bet1 = bet1_7
            bet2 = bet2_7
            s = s7
            F7 = fext * pow( pow(nu/nub,-s*bet1) + pow(nu/nub,-s*bet2), 1./s)
            ff = F7 * Ftilde8 * Ftilde9
   
        # Validity of spectrum 5
        if spectrum_case == 5:
            nuac = nu7
            fac = f7
            nusa = nu10
            fsa = f10
            nuc = nu11
            fc = f11
            num = nu9
            fm = f9
            # Prescription for the broadband soectra (formulae 4-9):
            Ftilde10 = pow(1. + pow(nu/nu10, s10 * (bet1_10-bet2_10)), -1./s10)
            Ftilde11 = pow(1. + pow(nu/nu11, s11 * (bet1_11-bet2_11)), -1./s11)
            Ftilde9 = pow(1. + pow(nu/nu9, s9 * (bet1_9 -bet2_9)), -1./s9 )
            nub = nu7
            fext = f7
            bet1 = bet1_7
            bet2 = bet2_7
            s = s7
            F7 = fext * pow( pow(nu/nub, -s*bet1) + pow(nu/nub, -s*bet2), -1./s)
            ff = F7 * Ftilde10 * Ftilde11 * Ftilde9
   
        #results = [ff, utils.mJy2mag(f, 'R'),spectrum_case,nuac,nusa,nuc,num,pls]
        #results = [ff,spectrum_case,nuac,nusa,nuc,num,pls]
        return ff

    @cython.boundscheck(False)
    def light_curve(self, double[:] time_series, double[:] frequencies):
        """ SImulate the afterglow light curve for a given wavelength"""
        #Make sure that we are working with arrays
        time_series = np.atleast_1d(time_series)
        frequencies = np.atleast_1d(frequencies)
        
        cdef int i, j
        cdef int t1 = len(time_series)
        cdef int nu1 = len(frequencies)
        cdef double[:,:] lc = np.zeros((t1,nu1)) 
        
        self.DL = cosmo.luminosity_distance(self.z).value
        
        # Start with wavelength as it is the biggest array most of the time
        for i in prange(nu1, nogil=True, num_threads=self.num_threads):
            for j in range(t1):
                # Flux in obs frame in mJy
                lc[j,i] = self.granot_sari(td=time_series[j],nu=frequencies[i])

        return np.array(lc)
