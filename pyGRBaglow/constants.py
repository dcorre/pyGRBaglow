"""
    Various constants used

    Unit abreviations are appended to the name, but powers are not
    specified. For instance, the gravitational constant has units "Mpc^3
    msun^-1 s^-2", but is called "G_const_Mpc_Msun_s".

    Most of these numbers are from google calculator.

"""
### If you add a constant, make sure to add a description too. ###

#Physical contants

H_planck = 6.6262e-34    # Planck constant (in J.s)
c_light_m_s = 2.9979e8   #299792458  # speed of ight in vacuum (in m/s)
K_boltzman = 1.3807e-23  # Boltzman constant in (J/K) 
e_elec = 1.60217646e-19  # charge of the electron (in C)
m_p = 1.67262158e-27     # mass of a proton (in kg)
m_e = 9.11e-31           # mass of an electron (in kg)
m_sun = 1.98892e30       # mass of the sun (in kg)

h_c =  H_planck * c_light_m_s # h*c
zero_celsius = 273.15    # Zero Celcius (K)

# Black body constants
BBT = H_planck / K_boltzman
BBO = 2 * H_planck / c_light_m_s**2.


# Frequencies
nu_XRT03 = 1e17     # Hz at 0.3 keV
nu_XRT3 = 1e18      # Hz at 3 keV
nu_XRT10 = 2e18     # Hz at 10 keV
nu_band_U = 8.3275e14                    # Hz
nu_band_B = 6.8134090909090912e+014      # Hz
nu_band_V = 5.4507272727272725e+014      # Hz
nu_band_R = 4.6842187500000000e+014      # Hz
nu_band_I = 3.7948101265822787e+014      # Hz
nu_band_J = 2.3792857142857144e+014      # Hz 
nu_band_H = 1.8736875000000000e+014      # Hz
nu_band_K = 1.3504054054054055e+014      # Hz


# Zero point magnitude
m0Jy_band_U = 1810;  # Jy mag=0
m0Jy_band_B = 4260;
m0Jy_band_V = 3640;
m0Jy_band_R = 3080;
m0Jy_band_I = 2550;
m0Jy_band_J = 1600;
m0Jy_band_H = 1080;
m0Jy_band_K = 670;

# --- magnitude AB offsets
""" AB magnitude system:
       This magnitude system is defined such that, 
       when monochromatic flux f is measured 
       in erg sec^-1 cm^-2 Hz^-1,
       with 1 Jy = 10^-23 erg sec^-1 cm^-2 Hz^-1
       m(AB) = -2.5 log(f) - 48.60
       offset is V = V(AB) + offAB_bandV
"""

offAB_band_U = -0.; # TO DO...
offAB_band_B = 0.163;
offAB_band_V = 0.044;
offAB_band_R = -0.055;
offAB_band_I = -0.309;
offAB_band_J = -0.; # TO DO...
offAB_band_H = -0.;
offAB_band_K = -0.;

# --- extinction laws
EXTINCTION_LAW_MW = 1
EXTINCTION_LAW_LMC = 2
#EXTINCTION_LAW_DOR30 = 3;
EXTINCTION_LAW_SMC = 3
EXTINCTION_LAW_LINEAR = 4
EXTINCTION_LAW_CALZETTI = 5
#EXTINCTION_LAW_FLAT1 = 5;
EXTINCTION_LAW_GRB1 = 6
EXTINCTION_LAW_GRB2 = 7


# Others
pc_cm = 3.085677e18  #3.08568025e18 # Parsec in cm
Mpc_cm = pc_cm * 1e6 # Megaparsec in cm
Mpc_km = Mpc_cm * 1e-5 #Megaparsec in km
angstrom_cm = 1e-8 # Angstrom in cm
yr_s = 365. * 24. * 60. * 60. # a year in s
Myr_s = 1.e6 * yr_s  # Megayear in s
Gyr_s = 1.e9 * yr_s # Gigayear in s
amu_g = 1.66053886e-24 # atomic mass unit in g
m_p_g = m_p * 1e-6 # mass of a proton in g
m_H_g = 1.00794 * amu_g # mass of a hydrogen atom in g
m_He_g = 4.002602 * amu_g #  mass of a helium atom in g
M_sun_g = 1.98892e33 # Solar mass in grams
c_light_km_s = c_light_m_s * 1e-3  # Speed of light in m/s
c_light_cm_s = c_light_m_s * 1e2 # Speed of light in cm/s
c_light_Mpc_s = c_light_cm_s / Mpc_cm # Speed of light in Mpc/s
c_light_Mpc_Gyr =  Gyr_s * c_light_cm_s / Mpc_cm  # Speed of light in Mpc/Gyr
H100_km_s_Mpc = 100. # km/s/Mpc
H100_s = 100. / Mpc_km # 100 km s^-1 Mpc^-1 in s^-1
G_const_Mpc_Msun_s = M_sun_g * (6.673e-8) / Mpc_cm**3. #Gravitational constant in Mpc^3 msun^-1 s^-2
lambda_Lya_0 = 1215.67 # Central wavelength of H Lyman-alpha in Angstroms
lambda_NV_0 = 1240.81 #Central wavelength of an NV doublet in Angstroms
alpha_B_cm_s_1e4 = 2.59e-13 #  hydrogen recombination coefficient at T=10^4 K
sigma_T_cm = 6.6524586e-25 # Thomson cross section in cm^2
sigma_T_Mpc = sigma_T_cm / (Mpc_cm ** 2.) # Thomson cross section in Mpc^2

