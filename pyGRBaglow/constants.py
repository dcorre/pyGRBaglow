"""
    Various constants used

    Unit abreviations are appended to the name, but powers are not
    specified. For instance, the gravitational constant has units "Mpc^3
    msun^-1 s^-2", but is called "G_const_Mpc_Msun_s".

    Most of these numbers are from google calculator.

    SHOULD MOVE TO ASTROPY.UNITS AT SOME POINTS!!!!!

"""
# If you add a constant, make sure to add a description too.

# Physical contants

# Planck constant (in J.s)
H_planck = 6.6262e-34
# speed of ight in vacuum (in m/s)
c_light_m_s = 2.9979e8
# Boltzman constant in (J/K)
K_boltzman = 1.3807e-23
# charge of the electron (in C)
e_elec = 1.60217646e-19
# mass of a proton (in kg)
m_p = 1.67262158e-27
# mass of an electron (in kg)
m_e = 9.11e-31
# mass of the sun (in kg)
m_sun = 1.98892e30

# h*c
h_c = H_planck * c_light_m_s
# Zero Celcius (K)
zero_celsius = 273.15

# Black body constants
BBT = H_planck / K_boltzman
BBO = 2. * H_planck / c_light_m_s**2.

# Zero point magnitude
# Jy mag=0
m0Jy_band_U = 1810
m0Jy_band_B = 4260
m0Jy_band_V = 3640
m0Jy_band_R = 3080
m0Jy_band_I = 2550
m0Jy_band_J = 1600
m0Jy_band_H = 1080
m0Jy_band_K = 670

# --- extinction laws
EXTINCTION_LAW_MW = 1
EXTINCTION_LAW_LMC = 2
EXTINCTION_LAW_SMC = 3
EXTINCTION_LAW_LINEAR = 4
EXTINCTION_LAW_CALZETTI = 5
EXTINCTION_LAW_GRB1 = 6
EXTINCTION_LAW_GRB2 = 7

# Others
# Parsec in cm
pc_cm = 3.085677e18
# Megaparsec in cm
Mpc_cm = pc_cm * 1e6
# Megaparsec in km
Mpc_km = Mpc_cm * 1e-5
# Angstrom in cm
angstrom_cm = 1e-8
# a year in s
yr_s = 365. * 24. * 60. * 60.
# Megayear in s
Myr_s = 1.e6 * yr_s
# Gigayear in s
Gyr_s = 1.e9 * yr_s
# atomic mass unit in g
amu_g = 1.66053886e-24
# mass of a proton in g
m_p_g = m_p * 1e-6
# mass of a hydrogen atom in g
m_H_g = 1.00794 * amu_g
#  mass of a helium atom in g
m_He_g = 4.002602 * amu_g
# Solar mass in grams
M_sun_g = 1.98892e33
# Speed of light in m/s
c_light_km_s = c_light_m_s * 1e-3
# Speed of light in cm/s
c_light_cm_s = c_light_m_s * 1e2
# Speed of light in Mpc/s
c_light_Mpc_s = c_light_cm_s / Mpc_cm
# Speed of light in Mpc/Gyr
c_light_Mpc_Gyr = Gyr_s * c_light_cm_s / Mpc_cm
# km/s/Mpc
H100_km_s_Mpc = 100.
# 100 km s^-1 Mpc^-1 in s^-1
H100_s = 100. / Mpc_km
# Gravitational constant in Mpc^3 msun^-1 s^-2
G_const_Mpc_Msun_s = M_sun_g * (6.673e-8) / Mpc_cm**3.
# Central wavelength of H Lyman-alpha in Angstroms
lambda_Lya_0 = 1215.67
# Central wavelength of an NV doublet in Angstroms
lambda_NV_0 = 1240.81
# hydrogen recombination coefficient at T=10^4 K
alpha_B_cm_s_1e4 = 2.59e-13
# Thomson cross section in cm^2
sigma_T_cm = 6.6524586e-25
# Thomson cross section in Mpc^2
sigma_T_Mpc = sigma_T_cm / (Mpc_cm ** 2.)
