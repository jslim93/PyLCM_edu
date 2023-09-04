import numpy as np

rho_a = 1.225E-3 #g/cm3
rho_liq = 1000.0 # density of liquid water (kg/m3)
rv = 461.0 # specific gas constant water vapor (J/kg/K)
l_v = 2.501E6 # latent heat of vaporization (J/kg)
cp = 1005.0 # specific heat of air at constant pressure (J/kg/K)
g = 9.81 # acceleration by gravity (m/s2)
r_a = 287.0 # specific gas constant air (J/kg/K)
pi = 3.141592653589 # Pi
molecularviscosity = 1.461E-5 # (dynamic) viscosity of air (kg/m/s)
l_s = 2.834E6 # latent heat of sublimation (J/kg)
muelq = 1.717e-05 # (dynamic) viscosity of air (kg/m/s)
p0 = 1000.0E2 # reference pressure (Pa)
molecular_weight_air = 0.028 # molecular weight of air (kg/mol)
molecular_weight_water = 0.018 # molecular weight of water (kg/mol)
rho_ice = 916.7 # desity ice (kg/m3)
k_boltzmann = 1.3806485279E-23 # (in J/K) Boltzmann constant
N_avogadro = 6.02214085774E23 # (1/mol) Avogadro number
sigma_ice_vap =  100.0E-3 # (in J/m2) surface tension ice/vapor interface (Pruppacher and Klett, p. 341)
nu_s = 1.0E13 #, &  ! frequency of vibration of water vapor molecule adsorbed on solid substrate (Pruppacher and Klett, p. 299)
h_planck = 6.62607004E-34 # (J*s) Planck's constant
R = 8.314459848 # universal gas konstant (J/(mol K))
m_nuc = 4.0 / 3.0 * pi * rho_ice * ( 1.0E-6 )**3 # initial ice size ! halving or doubling seems to have no impact
v_w  = molecular_weight_water / ( N_avogadro * rho_ice ) # volume of water molecule in ice
alpha_T = 0.96 # gas kinetic correction of thermal conduction (thermal accommodation coefficient)
alpha_m = 0.5 # gas kinetic correction of water vapor diffusion (deposition coefficient)
N_l =  N_avogadro * rho_liq / molecular_weight_water, # number of water molecules per volume of water
m_w = molecular_weight_water / N_avogadro # mass of a water molecule
n_s = 1.0E19 # concentration of water molecules in contact with surface of ice germ/nucleus; values from 5.85E18 to 1.0E19 are found in the literature
sigma_stefan_boltzmann = 5.670374419E-8 # (in W/m2/T4) Stefan Boltzmann konstant 
beta_env = 1.0E-5
vanthoff_aero = 2.00

#sea-salt aerosol
rho_aero = 2170.0 # density (kg/m3) 2.17 g/cm3
molecular_weight_aero = 0.058443 # (kg/mol)   58.443 g/mol

#-------------------------------------------------------
#parameter to draw spectra. 
#-------------------------------------------------------

alpha_spec   = 1.0    # bin spacing: mass(n+1) = mass(n) * 2^alpha_spec
r_start_spec = 0.005E-6 #! smallest radius for spectra; 0.01E-6 to see haze, cloud droplets, and rain; 1.5625 to see cloud droplets and rain
x_start_spec =4.0 / 3.0 * np.pi * rho_liq * ( r_start_spec )**3
r_end_spec = 8000.0E-6
x_end_spec = 4.0 / 3.0 * np.pi * rho_liq * ( r_end_spec )**3

n_bins  = int( np.log10( x_end_spec / x_start_spec ) / np.log10( 2.0**alpha_spec ) )

n_bins_spec = np.arange(1,n_bins)
xl = x_start_spec * (2**alpha_spec)**(n_bins_spec-1) # where is x_start_spec defined? Maybe just 1?
xr = x_start_spec * (2**alpha_spec)**(n_bins_spec)
xm = np.sqrt(xl * xr)
rm_spec = ( xm / ( 4.0 / 3.0 * np.pi * rho_liq ) )**(1.0/3.0)  # this is the mean radius of the bin
rl_spec = ( xl / ( 4.0 / 3.0 * np.pi * rho_liq ) )**(1.0/3.0)
rr_spec = ( xr / ( 4.0 / 3.0 * np.pi * rho_liq ) )**(1.0/3.0)

#-------------------------------------------------------
#parameter to microphysics. 
#-------------------------------------------------------

activation_radius_ts = 1.0E-6  # Activation radius in meters, example value
seperation_radius_ts = 25.0E-6