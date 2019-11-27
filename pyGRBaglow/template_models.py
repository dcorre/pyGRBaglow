# -*- coding: utf-8 -*-

"""Main module."""

import numpy as np
from collections import OrderedDict
from astropy.cosmology import Planck15 as cosmo
from pyGRBaglow import constants as cc
import gc
import sys


class Templates(object):
   """ Power law, broken power law..."""

   def __init__ (self, F0=500, t0=300,wvl0=6400):
       """
       F0: flux normaliqation in Jy
       t0: time corresponding to the normalisation, in s
       wvl0: wavelength corresponding to the normalisation, in angstrom
       """
       self.F0=F0
       self.t0=t0
       self.wvl0=wvl0

       return None

   def SPL(self,wvl,t,alpha,beta):
       """ Simple Power Law. wvl and t must bhave same dimensions as wvl0 and t0"""

       F= self.F0 * (t/self.t0)**(-alpha) * (wvl/self.wvl0)**(beta)
       return F

   def BPL(self,wvl,t,alpha1,alpha2,beta,s):
       """ Broken Power Law. wvl and t must bhave same dimensions as wvl0 and t0 """
       #print (self.F0,self.wvl0,self.t0)
       #print (wvl,t)
       #print (alpha1,alpha2,beta,s)
       F= self.F0 * (wvl/self.wvl0)**(beta) * ((t/self.t0)**(-s*alpha1) + (t/self.t0)**(-s*alpha2)) **(-1/s)
       return F

   def light_curve(self,wavelength,time,params,model='SPL'):
       """ build light cirves """
       time_series=np.atleast_1d(time)
       wavelength=np.atleast_1d(wavelength)
       lc=[]
       t1=len(time)
       wvl1=len(wavelength)
       for t in range(t1):
           sed=[]
           for wvl in range(wvl1):
               if model == 'SPL':
                   alpha,beta= params
                   SED=self.SPL(wavelength[wvl],time[t],alpha,beta)

               elif model == 'BPL':
                   alpha1,alpha2,beta,s=params
                   SED=self.BPL(wavelength[wvl],time[t],alpha1,alpha2,beta,s)

               sed.append( SED) # Flux in obs frame in same unit as F0

           lc.append(sed) # Flux in obs frame in same unit as F0
       return np.array(lc)
 
