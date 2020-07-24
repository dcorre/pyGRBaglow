from pyGRBaglow import constants as cc
import numpy as np

class reddening(object):
    """
    Class to to model the host and galactic dust and gas attenuation of a grb spectrum because 
    """

    def __init__(self,wavelength,z,Av):
         """
	 Parameters
         ----------
          Av : float (optional)
               Coefficient of attenuation (default = 0.1)
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
         NHI: hydrogen column density              
         """
         # So far we just assume that this is 4 four times higher than in our galaxy
         # Should be a free parameter if X ray data are used
         return 4.*self.Av * 1.79e21


    def extinction_laws_coeff(self,extinction_law):
         """ Define the parameters for the extinction law for each cases using the "Drude" model by Li and al 2008

         """
         if extinction_law.lower() == 'mw':
              c1=14.4
              c2=6.52
              c3=2.04
              c4=0.0519
         elif extinction_law.lower() == 'lmc':
              c1=4.47
              c2=2.39
              c3=-0.988
              c4=0.0221
         elif extinction_law.lower() == 'smc':
              c1=38.7
              c2=3.83
              c3=6.34
              c4=0.
         elif extinction_law.lower() == 'linear':
              c1=66.2
              c2=4.97
              c3=22.1
              c4=0.
         elif extinction_law.lower() == 'calzetti':
              c1=44.9
              c2=7.56
              c3=61.2
              c4=0.
         elif extinction_law.lower() == 'grb1':
              c1=0.025
              c2=0.048
              c3=-2.
              c4=0.
         elif extinction_law.lower() == 'grb2':
              c1=0.015
              c2=0.15
              c3=-2.
              c4=0.
         return [c1,c2,c3,c4]

    def Li07(self,extinction_law='smc',Xcut=False):
         """ Compute the optical depth due to dust reddening 

         Parameters
         ----------
         extinction_law : string
                          model of extinction to be used
                               - Milky Way
                               - LMC
                               - SMC (default)
                               - Linear
                               - Calzetti

         wavelength : float or array
                      wavelength in Angstrom

         Return 
         ------
         Alambda_over_Av: float or array
                          reddening law

         Trans_dust : float or array
                      Transmission coefficient due to dust reddening to apply to the spectrum
   
         """
         #Awavel_over_Av = np.ones(len(nus))
         #print (extinction_law)
         #Load the coefficient corresponding to the desired extinction curve template
         [c1, c2, c3, c4] = self.extinction_laws_coeff(extinction_law)

         nus=cc.c_light_m_s/(self.wavelength*1e-10)

         nu=nus*(1+self.z)
         wavel_mic = self.wavelength *1e-4 / (1+self.z)

         # Check if the wavelength is a scalar or an array
         # Otherwise can not iterate over 0-d array
         if not hasattr(wavel_mic, "__len__"): wavel_mic=np.array([wavel_mic])
        
         # First term accounting for the far UV extinction rise
         UV_extinction = c1 / ( (wavel_mic/0.08)**c2 + (0.08/wavel_mic)**c2 + c3 )

         # Second term accounting for the near-IR/visible extinction
         IR_vis_ext = ( 233.*(1. - c1/(6.88**c2 + 0.145**c2 + c3) - c4/4.60) ) / ( (wavel_mic/0.046)**2. + (0.046/wavel_mic)**2. + 90. )

         # Third term accounting for the 2175 Angstrom extinction bump
         PAH_bump = c4 / ( (wavel_mic/0.2175)**2. + (0.2175/wavel_mic)**2. - 1.95 )

         Awavel_over_Av = UV_extinction + IR_vis_ext + PAH_bump    # In the rest frame

         #Set arbitrarily the negative extinction to zero
         w=np.where(Awavel_over_Av < 0)
         Awavel_over_Av[w] = 0

         #Applied a cut for wavelength below 700 angstrom
         #Useful when coupling with Xray data
         if Xcut:
             w=np.where(wavel_mic < 0.07)
             Awavel_over_Av[w]=0 

         # Return optical depth due to dust reddening in funtion of wavelength
         Tau_dust = self.Av/1.086 * Awavel_over_Av

         Trans_dust = np.exp(-Tau_dust)

         w=np.where(Trans_dust<=0)
         Trans_dust[w]=0
         w=np.where(Trans_dust>1)
         Trans_dust[w]=1

         return [Awavel_over_Av,Trans_dust]
         

    def Pei92(self,Rv='none',law='smc',Xcut=False):
       wvl=self.wavelength*1e-4/(1+self.z)
       if law.lower() == 'smc' :
           if Rv == 'none': Rv=2.93
           a_1=185;a_2=27;a_3=0.005;a_4=0.010;a_5=0.012;a_6=0.03;
           wvl_1=0.042;wvl_2=0.08;wvl_3=0.22;wvl_4=9.7;wvl_5=18;wvl_6=25;
           b_1=90;b_2=5.50;b_3=-1.95;b_4=-1.95;b_5=-1.80;b_6=0.0;
           n_1=2.0;n_2=4.0;n_3=2.0;n_4=2.0;n_5=2.0;n_6=2.0;
 
           sum1=a_1/( (wvl/wvl_1)**n_1 + (wvl_1/wvl)**n_1 + b_1 )
           sum2=a_2/( (wvl/wvl_2)**n_2 + (wvl_2/wvl)**n_2 + b_2 )
           sum3=a_3/( (wvl/wvl_3)**n_3 + (wvl_3/wvl)**n_3 + b_3 )
           sum4=a_4/( (wvl/wvl_4)**n_4 + (wvl_4/wvl)**n_4 + b_4 )
           sum5=a_5/( (wvl/wvl_5)**n_5 + (wvl_5/wvl)**n_5 + b_5 )
           sum6=a_6/( (wvl/wvl_6)**n_6 + (wvl_6/wvl)**n_6 + b_6 )

       elif law.lower() == 'lmc':
           if Rv == 'none': Rv=3.16
           a_1=175;a_2=19;a_3=0.023;a_4=0.005;a_5=0.006;a_6=0.02;
           wvl_1=0.046;wvl_2=0.08;wvl_3=0.22;wvl_4=9.7;wvl_5=18;wvl_6=25;
           b_1=90;b_2=5.5;b_3=-1.95;b_4=-1.95;b_5=-1.8;b_6=0.0;
           n_1=2.0;n_2=4.5;n_3=2.0;n_4=2.0;n_5=2.0;n_6=2.0;

           sum1=a_1/( (wvl/wvl_1)**n_1 + (wvl_1/wvl)**n_1 + b_1 )
           sum2=a_2/( (wvl/wvl_2)**n_2 + (wvl_2/wvl)**n_2 + b_2 )
           sum3=a_3/( (wvl/wvl_3)**n_3 + (wvl_3/wvl)**n_3 + b_3 )
           sum4=a_4/( (wvl/wvl_4)**n_4 + (wvl_4/wvl)**n_4 + b_4 )
           sum5=a_5/( (wvl/wvl_5)**n_5 + (wvl_5/wvl)**n_5 + b_5 )
           sum6=a_6/( (wvl/wvl_6)**n_6 + (wvl_6/wvl)**n_6 + b_6 )

       elif law.lower() == 'mw':
           if Rv == 'none': Rv=3.08
           a_1=165;a_2=14;a_3=0.045;a_4=0.002;a_5=0.002;a_6=0.012;
           wvl_1=0.046;wvl_2=0.08;wvl_3=0.22;wvl_4=9.7;wvl_5=18;wvl_6=25;
           b_1=90;b_2=4.0;b_3=-1.95;b_4=-1.95;b_5=-1.8;b_6=0.0;
           n_1=2.0;n_2=6.5;n_3=2.0;n_4=2.0;n_5=2.0;n_6=2.0;

           sum1=a_1/( (wvl/wvl_1)**n_1 + (wvl_1/wvl)**n_1 + b_1 )
           sum2=a_2/( (wvl/wvl_2)**n_2 + (wvl_2/wvl)**n_2 + b_2 )
           sum3=a_3/( (wvl/wvl_3)**n_3 + (wvl_3/wvl)**n_3 + b_3 )
           sum4=a_4/( (wvl/wvl_4)**n_4 + (wvl_4/wvl)**n_4 + b_4 )
           sum5=a_5/( (wvl/wvl_5)**n_5 + (wvl_5/wvl)**n_5 + b_5 )
           sum6=a_6/( (wvl/wvl_6)**n_6 + (wvl_6/wvl)**n_6 + b_6 )

       #Need to check whether extrapolation is needed outside the range defined in Pei92
       Alambda_over_Ab= (1/Rv+1) * (sum1+sum2+sum3+sum4+sum5+sum6)
 
       #Applied a cut for wavelength below 700 angstrom
       #Useful when coupling with Xray data
       if Xcut:
           w=np.where(wvl < 0.07)
           Alambda_over_Ab[w]=0


       # Return optical depth due to dust reddening in funtion of wavelength
       Tau_dust = self.Av * Alambda_over_Ab / 1.086

       Trans_dust = np.exp(-Tau_dust)

       Trans_dust[Trans_dust < 0]=0
       Trans_dust[Trans_dust > 1]=1

       return [Alambda_over_Ab, Trans_dust]

    def gas_absorption(self, NHx=0.2):
         """ Compute the optical depth due to gas absorption

         Parameters
         ----------
         wavelength : array
               wavelength in Angstrom
        
         NHx: float
              Metal column density from soft Xrays absortpion (in units 1e22 /cm2/mag),
              expressed in units of equivalent hydrogen column density assuming solar abundances
              In Milky Way: NHx/Av = 1.7 to 2.2 1e21 cm-2/mag (default set to 2)

         Return
         ------

         Trans_gas : array
                   Transmission coefficient due to the gas absorption either in the galactic or host

         """

         nus = cc.c_light_m_s / (self.wavelength *1e-10)
         Tau_gas = np.zeros(len(nus))

         #Set NHx in units 1e22 /cm2/mag
         NHx*=1e22

         # Define the equivalent HI column density depending whether we study the host or our galaxy
         """
         if self.z == 0:
              NHI = self.Av * 1.79e21    # cm-2
         else:
              NHI = self.NHI_host()    # cm-2
         """
         NHI = self.Av * NHx

         for i in range(len(nus)):
              nu = nus[i] * (1+self.z)                   # photon frequency (Hz) in the rest frame
              E_kev = nu * cc.H_planck / (1e3 * cc.e_elec)    # photon energy (keV) in the rest frame
              E_kev2 = E_kev**2.
              E_kev3 = E_kev**3.

              #if E_kev < 13.6e-3:    #912 A (Lyman limit)
              #     c0=0; c1=0; c2=0
              if E_kev < 0.030:       # 41nm / 410A
                   c0=0.; c1=0.; c2=0.; edge='H'
              elif E_kev < 0.100:     # 12.4 nm
                   c0=17.3; c1=608.1; c2=-2150; edge='He';
              elif E_kev < 0.284:     #4.37 nm
                   c0=34.6; c1=267.9; c2=-476.1; edge='C';
              elif E_kev < 0.400:
                   c0=78.1; c1=18.8; c2=4.3; edge='N';
              elif E_kev < 0.532:
                   c0=71.4; c1=66.8; c2=-51.4; edge='O';
              elif E_kev < 0.707:
                   c0=95.5; c1=145.8; c2=-61.1; edge='Fe-L';
              elif E_kev < 0.867:
                   c0=308.9; c1=-380.6; c2=294.0; edge='Ne';
              elif E_kev < 1.303:
                   c0=120.6; c1=169.3; c2=-47.7; edge='Mg';
              elif E_kev < 1.840:
                   c0=141.3; c1=146.8; c2=-31.5; edge='Si';
              elif E_kev < 2.471:
                   c0=202.7; c1=104.7; c2=-17.0; edge='S';
              elif E_kev < 3.210:
                   c0=342.7; c1=18.7; c2=0.0; edge='Ar';
              elif E_kev < 4.038:
                   c0=352.2; c1=18.7; c2=0.0; edge='Ca';
              elif E_kev < 7.111:
                   c0=433.9; c1=-2.4; c2=0.75; edge='Fe';
              elif E_kev < 8.331:
                   c0=629.0; c1=30.9; c2=0.0; edge='Ni';
              elif E_kev < 10.:    # 124pm/1.24A
                   c0=701.2; c1=25.2; c2=0.0; edge='...';
              else:
                   c0=0.; c1=0.; c2=0.

              sige3 = (c0+c1*E_kev+c2*E_kev2) # Figure of M&M
              sig = sige3/E_kev3*1e-24        # cross section per hydrogen atom /cm2
              
              Tau_gas[i] = sig * NHI

         Trans_gas = np.exp(-Tau_gas)

         Trans_gas[Trans_gas < 1e-5] = 0
         Trans_gas[Trans_gas > 1] = 1

         return Trans_gas
