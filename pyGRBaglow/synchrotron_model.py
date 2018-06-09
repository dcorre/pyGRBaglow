# -*- coding: utf-8 -*-

"""Main module."""

import numpy as np
from collections import OrderedDict
from astropy.cosmology import Planck15 as cosmo
from . import constants as cc
import gc
import sys


class fireball_afterglow(object):
   """  Class to to model a grb afterglow according to the fireball model """


   def __init__(self,n0=0.001,eps_b=0.001,eps_e=0.001,E_iso=2e52,eta=0.1,p=2.4,Y=0,Mdot_loss=1e-5,Vw=1000,z=1,td=1,nu=1,ism_type=0,disp=0):
        
       self.parameters=OrderedDict()
       """
       #update settings with defaults
       self.parameters.update(dict(n0=0.001,
                                     eps_b=0.001,
                                     eps_e=0.4,
                                     E_iso=2e52,
                                     eta=0.1,
                                     p=2.4,
                                     Y=0,
                                     Mdot_loss=1e-5,
                                     Vw=1000,
                                     z=1.,
                                     td=1.,
                                     nu=1.,
                                     ism_type=0,
                                     disp=0))

       """

       #Define parameters
       self.parameters['n0'] = float(n0)
       self.parameters['eps_b'] = float(eps_b)
       self.parameters['eps_e'] = float(eps_e)
       self.parameters['E_iso'] = float(E_iso)
       self.parameters['eta']= float(eta)
       self.parameters['p'] = float(p)
       self.parameters['Y'] = float(Y)
       self.parameters['Mdot_loss'] = float(Mdot_loss)
       self.parameters['Vw'] = float(Vw)
       self.parameters['z']=float(z)
       self.parameters['td']=float(td)
       self.parameters['nu']=float(nu)
       self.parameters['ism_type']=int(ism_type)
       self.parameters['disp']=int(disp)

   def A(self):
       """
       No idea

       Parameters
       ----------
       none              

       """
       return (self.parameters['Mdot_loss'] * cc.m_sun * 1e3 /365.25/86400.) / (4*np.pi*self.parameters['Vw']*1e5)


   def granot_sari(self,td=1,nu=1e15):
       """ From Granot et Sari """
      
       z  = self.parameters['z']
       #print (hex(id(z)))
       #can not call astropy here. For some reason it uses a lot of RAM. When called millions of times it is using too much GB. So it is called outside. Might be solved by updating astropy version.
       #DL = cosmo.luminosity_distance(self.parameters['z']).value
       DL=self.DL
       #DL=25
       #td = self.parameters['td']
       #nu = self.parameters['nu']
       ism_type = self.parameters['ism_type']
       disp = self.parameters['disp']

       t_obs = td * 86400.  # duration in seconds

       #Load the parameters
       eps_b = self.parameters['eps_b']
       eps_e = self.parameters['eps_e']
       E52 = (1-self.parameters['eta'])/self.parameters['eta'] * self.parameters['E_iso']*1e-52
       #d = DL*cc.Mpc_cm / (1e28)
       # Divide by 2e28 to be consistent with Turpin's catalogue.
       d = DL*cc.Mpc_cm / (1e28)
       #d = 4.9190e+28 / (2e28)
       n0 = self.parameters['n0']
       p = self.parameters['p']
       eps_ebar = eps_e *(p-2)/(p-1)
       Y = self.parameters['Y']

       nu14 = nu /1e14

       # Constant ISM
       if ism_type == 0:
   
            if disp == 1:
                 print ("CONSTANT ISM: \n")
                 print ('Forward shock, at tobs='+str(t_obs)+' seconds\n')
   
            # Granot et Sari 2002 (Section 5, page 825)
            evolution_param = n0 * E52**(4./7) * eps_b**(9./7)
            if evolution_param < 18:
                 evolution_case = 512
            else:
                 evolution_case = 432
   
            # Granot et Sari 2002 (Table 2)
            nu1=1.24*(p-1)**(3./5)/(3.*p+2)**(3./5)*1e9*(1+z)**(-1.)*eps_ebar**(-1)*eps_b**(1./5)*n0**(3./5)*E52**(1./5)
            nu2=3.73*(p-0.67)*1e15*(1+z)**(1./2)*E52**(1./2)*eps_ebar**(2.)*eps_b**(1./2)*td**(-3./2)
            nu3=6.37*(p-0.46)*1e13*np.exp(-1.16*p)*(1+z)**(-1./2)*eps_b**(-3./2)*n0**(-1.)*E52**(-1./2)*td**(-1./2)
            nu4=5.04*(p-1.22)*1e16*(1+z)**(1./2)*eps_ebar**(2.)*eps_b**(1./2)*E52**(1./2)*td**(-3./2)
            nu5=3.59*(4.03-p)*1e9*np.exp(2.34*p)*( eps_ebar**(4*(p-1))*eps_b**(p+2)*n0**(4.)*E52**(p+2)/(1+z)**(6-p)/td**(3*p+2) )**(1./(2*(p+4)))
            nu6=3.23*(p-1.76)*1e12*( eps_ebar**(4.*(p-1))*eps_b**(p-1)*n0**(2.)*E52**(p+1)/(1+z)**(7.-p)/td**(3.*(p+1)))**(1./(2*(p+5)))
            nu7=1.12*(3*p-1)**(8/5)/(3*p+2)**(8./5)*1e8*(1+z)**(-13./10)*eps_ebar**(-8./5)*eps_b**(-2./5)*n0**(3./10)*E52**(-1./10)*td**(3./10)
            nu8=1.98*1e11*(1+z)**(-1./2)*n0**(1./6)*E52**(1./6)*td**(-1./2)
            nu9=3.94*(p-0.74)*1e15*(1+z)**(1./2)*eps_ebar**(2.)*eps_b**(1./2)*E52**(1./2)*td**(-3./2)
            nu10=1.32*1e10*(1+z)**(-1./2)*eps_b**(6./5)*n0**(11./10)*E52**(7./10)*td**(-1./2)
            nu11=5.86*1e12*(1+z)**(-1./2)*eps_b**(-3./2)*n0**(-1.)*E52**(-1./2)*td**(-1./2)
   
            f1=0.647*(p-1)**(6./5)/( (3.*p-1)*(3.*p+2)**(1./5) ) * (1+z)**(1./2)*eps_ebar**(-1.)*eps_b**(2./5)*n0**(7./10)*E52**(9./10)*td**(1./2)*d**(-2.)
            f2=9.93*(p+0.14)*(1+z)*eps_b**(1./2)*n0**(1./2)*E52*d**(-2.)
            f3=4.68*np.exp(4.82*(p-2.5))*1e3*(1+z)**((p+1)/2.)*eps_ebar**(p-1)*eps_b**(p-1./2)*n0**(p/2.)*E52**((p+1)/2.)*td**((1-p)/2.)*d**(-2.)
            f4=3.72*(p-1.79)*1e15*(1+z)**(7./2)*eps_ebar**(5.)*eps_b*n0**(-1./2)*E52**(3./2)*td**(-5./2)*d**(-2.)
            f5=20.8*(p-1.53)*np.exp(2.56*p)*d**(-2.)*( (1+z)**(7.*p+3)*eps_b**(2.*p+3)*E52**(3.*p+7)/eps_ebar**(10.*(1-p))/td**(5*(p-1)) )**(1./(2.*(p+4)))
            f6=76.9*(p-1.08)*np.exp(2.06*p)*d**(-2.)*( (1+z)**(7.*p+5)*eps_b**(2.*p-5)*E52**(3.*p+5)/eps_ebar**(10.*(1-p))/n0**(p)/td**(5.*(p-1)) )**(1./(2.*(p+5)))
            f7=5.27*(3*p-1)**(11./5)/( (3.*p+2)**(11./5) )*1e-3 * (1+z)**(-1./10)*eps_ebar**(-11./5)*eps_b**(-4./5)*n0**(1./10)*E52**(3./10)*td**(11./10)*d**(-2.)
            f8=154*(1+z)*eps_b**(-1./4)*n0**(-1./12)*E52**(2./3)*d**(-2.)
            f9=0.221*(6.27-p)*(1+z)**(1./2)*eps_ebar**(-1.)*eps_b**(-1./2)*E52**(1./2)*td**(1./2)*d**(-2.)
            f10=3.72*(1+z)*eps_b**(7./5)*n0**(6./5)*E52**(7./5)*d**(-2.)
            f11=28.4*(1+z)*eps_b**(1./2)*n0**(1./2)*E52*d**(-2.)
   
   
            #print (hex(id(eps_ebar)))
            s1=1.64
            s2=1.84-0.40*p
            s3=1.15-0.06*p
            s4=3.44*p-1.41
            s5=1.47-0.21*p
            s6=0.94-0.14*p
            s7=1.99-0.04*p
            s8=0.907
            s9=3.34-0.82*p
            s10=1.213
            s11=0.597
   
            # Compute the choise of the spectrum case
            spectrum_case=0
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
            if disp:
                 print ("WINDY environment: \n")
                 print ('Forward shock, at tobs='+str(t_obs)+ ' seconds\n')
   
            A_star = self.A() /(5e11)
   
            evolution_param = A_star * eps_ebar**(-1.) * E52**(-3./7) * eps_b**(2./7)
   
            if evolution_param > 100:
                 evolution_case = 4512
            else:
                 evolution_case = 432
   
   
            # Granot et Sari 2002 (Table 2)
            nu1=8.31*(p-1)**(3./5)/(3.*p+2)**(3./5)*1e9*(1+z)**(-2./5)*eps_ebar**(-1.)*eps_b**(1./5)*A_star**(6./5)*E52**(-2./5)*td**(-3./5)
            nu2=4.02*(p-0.69)*1e15*(1+z)**(1./2)*E52**(1./2)*eps_ebar**(2.)*eps_b**(1./2)*td**(-3./2)
            nu3=4.40*(3.45-p)*1e10*np.exp(0.45*p)*(1+z)**(-3./2)*eps_b**(-3./2)*A_star**(-2.)*E52**(1./2)*td**(1./2)
            nu4=8.08*(p-1.22)*1e16*(1+z)**(1./2)*eps_ebar**(2.)*eps_b**(1./2)*E52**(1./2)*td**(-3./2)
            nu5=1.58*(4.10-p)*1e10*np.exp(2.16*p)*( eps_ebar**(4.*(p-1))*eps_b**(p+2)*A_star**(8.)/(1+z)**(2-p)/E52**(2-p)/td**(3.*(p+2)) )**(1./(2.*(p+4)))
            nu6=4.51*(p-1.73)*1e12*( eps_ebar**(4.*(p-1))*eps_b**(p-1)*A_star**(4.)*E52**(p-1)/(1+z)**(5-p)/td**(3.*p+5) )**(1./(2.*(p+5)))
            nu7=1.68*(3*p-1)**(8./5)/(3.*p+2)**(8./5)*1e8*(1+z)**(-1.)*eps_ebar**(-8./5)*eps_b**(-2./5)*A_star**(3./5)*E52**(-2./5)
            nu8=3.15*1e11*(1+z)**(-1./3)*A_star**(1./3)*td**(-2./3)
            nu9=3.52*(p-0.31)*1e15*(1+z)**(1./2)*eps_ebar**(2.)*eps_b**(1./2)*E52**(1./2)*td**(-3./2)
            # Following line must be updated
            nu10=1.32*1e10*(1+z)**(-1./2)*eps_b**(6./5)*n0**(11./10)*E52**(7./10)*td**(-1./2)
            # Following line must be updated
            nu11=5.86*1e12*(1+z)**(-1./2)*eps_b**(-3./2)*n0**(-1.)*E52**(-1./2)*td**(-1./2)
   
            f1=9.19*(p-1)**(6./5)/( (3.*p-1)*(3*p+2)**(1./5) ) * (1+z)**(6./5)*eps_ebar**(-1.)*eps_b**(2./5)*A_star**(7./5)*E52**(1./5)*td**(-1./5)*d**(-2.)
            f2=76.9*(p+0.12)*(1+z)**(3./2)*eps_b**(1./2)*A_star*E52**(1./2)*td*(-1./2)*d**(-2.)
            f3=8.02*np.exp(7.02*(p-2.5))*1e5*(1+z)**(p+1./2)*eps_ebar**(p-1)*eps_b**(p-1./2)*A_star**(p)*E52**(1./2)*td**(1./2-p)*d**(-2.)
            f4=3.04*(p-1.79)*1e15*(1+z)**(3.)*eps_ebar**(5.)*eps_b*A_star**(-1.)*E52**(2.)*td**(-2.)*d**(-2.)
            f5=158*(p-1.48)*np.exp(2.24*p)*d**(-2.)*( (1+z)**(6.*p+9)*eps_b**(2.*p+3)*E52**(4.*p+1)/eps_ebar**(10.*(1-p))/A_star**(2.*(p-6))*td**(4*p+1) )**(1./(2*(p+4)))
            f6=78.6*(p-1.12)*np.exp(1.89*p)*d**(-2.)*( (1+z)**(6.*p+5)*eps_b**(2.*p-5)*E52**(4.*p+5)/eps_ebar**(10*(1-p))/A_star**(2.*p)/td**(4.*p-5) )**(1./(2.*(p+5)))
            f7=3.76*(3*p-1)**(11./5)/( (3.*p+2)**(11./5) )*1e-3 * eps_ebar**(-11./5)*eps_b**(-4./5)*A_star**(1./5)*E52**(1./5)*td*d**(-2.)
            f8=119*(1+z)**(11./12)*eps_b**(-1./4)*A_star**(-1./6)*E52**(3./4)*td**(1./12)*d**(-2.)
            f9=0.165*(7.14-p)*(1+z)**(1./2)*eps_ebar**(-1.)*eps_b**(-1./2)*E52**(1./2)*td**(1./2)*d**(-2.)
            # Following line must be updated
            f10=3.72*(1+z)*eps_b**(7./5)*n0**(6./5)*E52**(7./5)*d**(-2.)
            # Following line must be updated
            f11=28.4*(1+z)*eps_b**(1./2)*n0**(1./2)*E52*d**(-2.)
   
            s1=1.06
            s2=1.76-0.38*p
            s3=0.80-0.03*p
            s4=3.63*p-1.60
            s5=1.25-0.18*p
            s6=1.04-0.16*p
            s7=1.97-0.04*p
            s8=0.893
            s9=3.68-0.89*p
            # Following line must be updated
            s10=1.213
            # Following line must be updated
            s11=0.597
   
            #Compute the choise of the spectrum case
            spectrum_case=0
   
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
   
       if disp:
            print ('Evolution case : %s \n' % evolution_case)
            print ('Spectrum Type : %s \n' % spectrum_case)
   
       # Inverse Compton discussed in page 827 
       # Y = eps_e /eps_b if eps_e << eps_b
       # Y = sqrt(eps_e/eps_b) if eps_e >> eps_b (Sari, Piran et Narayan 1996)
   
       if Y > 0. and eps_b < 10.*eps_e:
            if disp:
                 print ('Inverse Compton occurs')
   
            nu1=nu1*(1+Y)
            nu2=nu2*(1+Y)
            nu3=nu3*(1+Y)
            nu4=nu4*(1+Y)
            nu5=nu5*(1+Y)
            nu6=nu6*(1+Y)
            nu7=nu7*(1+Y)
            nu8=nu8*(1+Y)
            nu9=nu9*(1+Y)
            nu10=nu10*(1+Y)
            nu11=nu11*(1+Y)
   
       else:
            if disp:
                 print ('No inverse Compton expected')
   
       # Granot et SAri 2002 (Table 2)
       bet1_1=2.
       bet1_2=1./3
       bet1_3=(1-p)/2.
       bet1_4=2.
       bet1_5=5./2
       bet1_6=5./2
       bet1_7=2.
       bet1_8=11./8
       bet1_9=-1./2
       bet1_10=11./8
       bet1_11=1./3
   
       bet2_1=1./3
       bet2_2=(1-p)/2.
       bet2_3=-p/2.
       bet2_4=5./2
       bet2_5=(1-p)/2.
       bet2_6=-p/2.
       bet2_7=11./8
       bet2_8=-1./2
       bet2_9=-p/2.
       bet2_10=1./3
       bet2_11=-1./2
   
       nuac=0.
       nusa=0.
       nuc=0.
       num=0.
       f=0.
       ff=0.
   
   
       # Validity of Spectrum 1
       if spectrum_case == 1:
            nusa=nu1
            fsa=f1
            num=nu2
            fm=f2
            nuc=nu3
            fc=f3
            # Prescription for the broad band spectra (formulae 4-9)
            Ftilde2 =(1+(nu/nu2 )**(s2 *(bet1_2 -bet2_2 )))**(-1./s2 )
            Ftilde3 =(1+(nu/nu3 )**(s3 *(bet1_3 -bet2_3 )))**(-1./s3 )
            nub=nu1
            fext=f1
            bet1=bet1_1
            bet2=bet2_1
            s=s1
            F1=fext*( (nu/nub)**(-s*bet1) + (nu/nub)**(-s*bet2) )**(-1./s)
            ff = F1*Ftilde2*Ftilde3
   
       # Validity of spectrum 2
       if spectrum_case == 2:
            num=nu4
            fm=f4
            nusa=nu5
            fsa=f5
            nuc=nu3
            fc=f3
            # Prescription for the broadband spectra (formulae 4-9)
            Ftilde5 =(1+(nu/nu5 )**(s5 *(bet1_5 -bet2_5 )))**(-1/s5 )
            Ftilde3 =(1+(nu/nu3 )**(s3 *(bet1_3 -bet2_3 )))**(-1/s3 )
            nub=nu4
            fext=f4
            bet1=bet1_4
            bet2=bet2_4
            s=s4
            if nub == nu4:
                 phi4=(nu/nu4) # (3)
                 F4=f4*(phi4*phi4*np.exp(-s*phi4**(2./3))+phi4**(5./2))
            else:
                 F4=fext*( (nu/nub)**(-s*bet1) + (nu/nub)**(-s*bet2) )**(-1./s)
            ff=F4*Ftilde5*Ftilde3
   
   
       # Validity of spectrum 3
       if spectrum_case == 3:
            num=nu4
            fm=f4
            nusa=nu6
            fsa=f6
            # Prescription for the broadband spectra (formulae 4-9)
            Ftilde6 =(1+(nu/nu6 )**(s6 *(bet1_6 -bet2_6 )))**(-1./s6 )
            nub=nu4
            fext=f4
            bet1=bet1_4
            bet2=bet2_4
            s=s4
            if nub == nu4:
                 phi4=(nu/nu4) # (3)
                 F4=f4*(phi4*phi4*np.exp(-s*phi4**(2./3))+phi4**(5./2))
            else:
                 F4=fext*( (nu/nub)**(-s*bet1) + (nu/nub)**(-s*bet2) )**(-1./s)
            ff=F4*Ftilde6
   
   
       # Validity of spetrum 4 
       if spectrum_case == 4:
            nuac=nu7
            fac=f7
            nusa=nu8
            fsa=f8
            num=nu9
            fm=f9
            # Prescription for the broadband spectra (formulae 4-9)
            Ftilde8 =(1+(nu/nu8 )**(s8 *(bet1_8 -bet2_8 )))**(-1./s8 )
            Ftilde9 =(1+(nu/nu9 )**(s9 *(bet1_9 -bet2_9 )))**(-1./s9 )
            nub=nu7
            fext=f7
            bet1=bet1_7
            bet2=bet2_7
            s=s7
            F7=fext*( (nu/nub)**(-s*bet1) + (nu/nub)**(-s*bet2) )**(-1./s)
            ff=F7*Ftilde8*Ftilde9
   
       # Validity of spectrum 5
       if spectrum_case == 5:
            nuac=nu7
            fac=f7
            nusa=nu10
            fsa=f10
            nuc=nu11
            fc=f11
            num=nu9
            fm=f9
            # Prescription for the broadband soectra (formulae 4-9):
            Ftilde10=(1+(nu/nu10)**(s10*(bet1_10-bet2_10)))**(-1./s10)
            Ftilde11=(1+(nu/nu11)**(s11*(bet1_11-bet2_11)))**(-1./s11)
            Ftilde9 =(1+(nu/nu9 )**(s9 *(bet1_9 -bet2_9 )))**(-1./s9 )
            nub=nu7
            fext=f7
            bet1=bet1_7
            bet2=bet2_7
            s=s7
            F7=fext*( (nu/nub)**(-s*bet1) + (nu/nub)**(-s*bet2) )**(-1./s)
            ff=F7*Ftilde10*Ftilde11*Ftilde9
   
       #results = [ff, utils.mJy2mag(f, 'R'),spectrum_case,nuac,nusa,nuc,num,pls]
       #results = [ff,spectrum_case,nuac,nusa,nuc,num,pls]
       return ff


   def light_curve(self,time_series,frequencies):
       """ SImulate the afterglow light curve for a given wavelength"""
       #Make sure that we are working with arrays
       time_series=np.atleast_1d(time_series)
       frequencies=np.atleast_1d(frequencies)


       self.DL = cosmo.luminosity_distance(self.parameters['z']).value

       lc=[] 
       t1=len(time_series)
       nu1=len(frequencies)
       for t in range(t1):
           # Check if nu is a sequence or float
           sed=[]
           for nu in range(nu1):
               sed.append( self.granot_sari(td=time_series[t],nu=frequencies[nu]) ) # Flux in obs frame in mJy 
  
           lc.append(sed) # Flux in obs frame in mJy

       return np.array(lc)
