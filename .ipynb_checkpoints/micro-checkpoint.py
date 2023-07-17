import numpy as np
from matplotlib import pyplot as plt
from parameters import *
import random

class particles:
    def __init__(self,n):
        self.id     = n
        self.M      = 1.0
        self.A      = 1.0 #weighting factor
        self.r      = 1.0
    def __shuffle__(self):
        random.shuffle(self.particles)
        
def r_equi(S,T,r_aerosol, rho_aero,molecular_weight_aero):
    S_internal = min( S, -0.0001 )
    afactor = 2.0 * sigma_air_liq(T) / ( rho_liq * rv * T ) # curvature effect
    bfactor = vanthoff_aero * rho_aero * molecular_weight_water / ( rho_liq * molecular_weight_aero ) # solute effect

#     Iterative solver ( 0 = S - A / r + B / r^3 => r = ( B / ( A / r - S ) )^(1/3) )
    r_equi_0 = 1.0
    r_equi   = 1.0E-6
    
    while ( abs( ( r_equi - r_equi_0 ) / r_equi_0 ) > 1.0E-20 ):
        r_equi_0 = r_equi
        r_equi   = ( ( bfactor * r_aerosol**3 ) / ( afactor / r_equi - S_internal ) )**(1.0/3.0)
    return(r_equi)

def sigma_air_liq(tabs):
    tabs_c = tabs - 273.15
    #!
    #!--    Pruppacher and Klett (1997), Eq. 5-12
    sigma_air_liq = 75.93 + 0.115 * tabs_c    + 6.818e-2 * tabs_c**2 + 6.511e-3 * tabs_c**3 + 2.933e-4 * tabs_c**4 + 6.283e-6 * tabs_c**5 + 5.285e-8 * tabs_c**6
    sigma_air_liq = sigma_air_liq * 1.0E-3
    
    return(sigma_air_liq)

def radius_liquid_euler_old(r_ini,dt_int,r0,G_pre,supersat,ventilation_effect,afactor,bfactor,r_aero,D_pre,radiation):   
    
    r_eul     = r_ini
    r_eul_old = r_ini

    dt_eul = dt_int
    t_eul  = 0.0

    while ( t_eul < dt_int - 1.0E-20 ):

        for m in range(500):

            dr2dt      = 2.0 * G_pre * ventilation_effect * ( supersat - afactor / r_eul + bfactor * r_aero**3 / r_eul**3 - D_pre * radiation * r_eul ) * r_eul / ( r_eul + r0 )
            d2r2dtdr2  = G_pre * ventilation_effect * ( afactor * r_eul**3 - bfactor * r_aero**3 * ( 3.0 * r_eul + 2.0 * r0 ) - r_eul**3 * ( D_pre * radiation * r_eul * ( r_eul + 2.0 * r0 ) - r0 * supersat ) ) / ( r_eul**4 * ( r_eul + r0 )**2 )

            dt_eul = min( 0.5 * abs( 1.0 / d2r2dtdr2 ), dt_int - t_eul, dt_int )

        #-- To speed up the Newton-Raphson scheme, the square root is executed at every iteration
            f       = r_eul**2 - r_ini**2 - dt_eul * dr2dt
            dfdr2   = 1.0 - dt_eul * d2r2dtdr2 #! = 1.0 - ( 0.0 + dt_eul * d2r2dtdr2 )

            r_eul   = np.sqrt( max( r_eul**2 - f / dfdr2, r_aero**2 ) ) # Newton-Raphson scheme (2nd order)

            rel_change = abs( r_eul - r_eul_old ) / r_eul_old
            r_eul_old  = r_eul

            if rel_change < 1.0E-12 :
                break

        t_eul = t_eul + dt_eul
        r_ini = r_eul

    return r_ini


import numpy as np

def radius_liquid_euler_py(r_ini, dt_int, r0, G_pre, supersat, ventilation_effect, afactor, bfactor, r_aero, D_pre, radiation):
    r_eul = r_ini
    r_eul_old = r_ini

    dt_eul = dt_int
    t_eul = 0.0

    while t_eul < dt_int - 1.0E-20:
        for m in range(500):
            r_eul_plus_r0 = r_eul + r0
            r_eul_squared = r_eul**2
            r_eul_cubed = r_eul**3
            r_aero_cubed = r_aero**3

            dr2dt = 2.0 * G_pre * ventilation_effect * (supersat - afactor / r_eul + bfactor * r_aero_cubed / r_eul_cubed - D_pre * radiation * r_eul) * r_eul / r_eul_plus_r0
            d2r2dtdr2 = G_pre * ventilation_effect * (afactor * r_eul_cubed - bfactor * r_aero_cubed * (3.0 * r_eul + 2.0 * r0) - r_eul_cubed * (D_pre * radiation * r_eul * (r_eul + 2.0 * r0) - r0 * supersat)) / (r_eul**4 * r_eul_plus_r0**2)

            f = r_eul_squared - r_ini**2 - dt_eul * dr2dt
            dfdr2 = 1.0 - dt_eul * d2r2dtdr2

            r_eul_squared = max(r_eul_squared - f / dfdr2, r_aero**2)
            r_eul = np.sqrt(r_eul_squared)

            rel_change = abs(r_eul - r_eul_old) / r_eul_old
            r_eul_old = r_eul

            if rel_change < 1.0E-12:
                break

        t_eul += dt_eul
        r_ini = r_eul

    return r_ini
