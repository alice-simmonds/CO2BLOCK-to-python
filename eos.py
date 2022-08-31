# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 16:55:10 2022

@author: simmonds
"""

# This function returns brine viscosity and CO2 density and viscosity.
# Brine viscosity is calculated according to Batzle and Wang (1992)
# CO2 density is calculated according to Redlich and Kwong (1949), with the parameters proposed by Spycher et al. (2003)
# CO2 viscosity is calculated acccording to Altunin and Sakhabeetdinov (1972)

import numpy as np
def eos(T0_mean,p,salinity,co2dens):                    # T in deg C, p in MPa, salinity in ppm/1e6
    
    #brine viscosity
    T = T0_mean
    brineviscosity = (0.1+0.333*salinity+(1.65+91.9*salinity**3)*np.exp(-(0.42*(salinity**0.8-0.17)**2+0.045)*T**0.8))/1e3
    T = T0_mean+273.15                      # convert C to K
    p = p*1e6                               # convert MPa to Pa
    
    # CO2 density
    a0 = 7.54                               # constant [Pa m6 K^0.5 mol^-2]
    a1 = -4.13e-3                           # constant [Pa m6 K^0.5 mol^-2]
    B = 2.78e-5                             # constant [m3/mol]
    a = a0+a1*T
    R = 8.314472                            # gas constant [m3 Pa K^-1 mol^-1]
    
    
    coeff1 = 1
    coeff2 = -(R*T/p)
    coeff3 = -(R*T*B/p-a/p/np.sqrt(T)+B**2)
    coeff4 = -(a*B/p/np.sqrt(T))
    coeff = [coeff1,coeff2,coeff3,coeff4]
    V = (np.roots(coeff))                   # molar volume [m3/mol]
    Vreal = np.real(V[np.isreal(V)])
    co2density = 0.044/Vreal                # [kg/m3]
    
    
    #CO2 viscosity
    a10 = 0.248566120
    a11 = 0.004894942
    a20 = -0.373300660
    a21 = 1.22753488
    a30 = 0.363854523
    a31 = -0.774229021
    a40 = -0.0639070755
    a41 = 0.142507049
    Tr = T/304                              # reduced temperature
    dens_r = co2dens/468                    # reduced density
    mu_0 = Tr**0.5 * (27.2246461 - 16.6346068 / Tr + 4.66920556 / (Tr**2)) * 1e-6
    co2viscosity = mu_0*np.exp(a10*dens_r + a11*dens_r/Tr +a20*dens_r**2 + a21*dens_r**2/Tr + a30*dens_r**3 + a31*dens_r**3/Tr + a40*dens_r**4 + a41*dens_r**4/Tr)         #[Pa/s]
    
    return [brineviscosity,co2density,co2viscosity]
