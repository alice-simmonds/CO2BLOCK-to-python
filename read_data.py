# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 14:43:39 2022

@author: simmonds
"""

import pandas as pd
import numpy as np


#read parameters and evaluate maximum sustainable pressure
def read_data(data):
 
    #default parameters (used in the case they aren't provided)
   litho_grad = 23                      # lithostatic gradient [MPa/km]
   hydro_grad = 10                      # hydrostatic gradient [MPa/km]
   temp_grad = 33                       # temperature gradient [C/km]
   def_k0 = 0.7                         # default stress ratio s3/s1
   def_friction_angle = 30              # default rock friction angle
   def_cohesion = 0                     # default rock cohesion [MPa]
   def_cr = 5*1e-4                      # default rock compressibility [MPa^-1]
   def_cw = 3*1e-4                      # deault water compressibility [MPa^-1]
   def_salinity = 180000                # default salinity [ppm]
    
    
    #read data in from Excel
    
   domain_type = data.values[0,1]       # domain confinement
   depth = data.values[0,2]             # shallowest depth of reservoir [m]
   depth_mean = data.values[0,3]        # mean depth of reservoir [m]
   thick = data.values[0,4]             # thickness of reservoir [m]
   area_res = data.values[0,5]          # area of reservoir [km^2]
   perm = data.values[0,6]*1e-15        # intrinsic permeability [m^2]
   por = data.values[0,7]               # porosity [-]
   cr = data.values[0,8]/1e6            # rock compressibility [1/Pa]
   cw = data.values[0,9]/1e6            # water compressibility [1/Pa]
   dens_c = data.values[0,10]*1e3       # density of CO2 [kg/m^3]
   visc_c = data.values[0,11]/1e3       # viscosity of CO2 [Pa.s]
   visc_w = data.values[0,12]/1e3       # viscosity of water [Pa.s]
   pres0 = data.values[0,13]            # pressure at the top of the reservoir [MPa]
   pres0_mean = data.values[0,14]       # pressure at the center of the reservoir [MPa]
   T0_mean = data.values[0,15]          # temperature at centre of reservoir [C]
   salinity = data.values[0,16]/1e6     # aquifer salinity [ppm/1e6]
   s1_tot = data.values[0,17]           # total maximum principal stress at the reservoir [MPa]
   stress_ratio = data.values[0,18]     # ratio of principal effective stresses (s3/s1) [-]
   friction = data.values[0,19]         # rock friction angle [deg]
   cohesion = data.values[0,20]         # rock cohesion coefficient [MPa]
   tens_strength = data.values[0,21]    # rock tensile stength [MPa]
   
   
   if (domain_type == "Open") or (domain_type == "open"):
       rc = np.inf
   if (domain_type == "Closed") or (domain_type == "closed"):
       rc = np.sqrt((area_res*1e6)/np.pi)


   from eos import eos       
#calculate some parameters if not given
   if pres0 == 0 or np.isnan(pres0):
       pres0 = hydro_grad*depth/1000
   if pres0_mean == 0 or np.isnan(pres0_mean):
       pres0_mean = hydro_grad*depth_mean/1000
   if T0_mean == 0 or np.isnan(T0_mean):
       T0_mean = temp_grad*depth_mean/1000 + 15
   if s1_tot == 0 or np.isnan(s1_tot):
       s1_tot = litho_grad*depth/1000
   if stress_ratio == 0 or np.isnan(stress_ratio):
       stress_ratio = def_k0
   if friction == 0 or np.isnan(friction):
       friction = def_friction_angle
   if cohesion == 0 or np.isnan(cohesion):
       cohesion = def_cohesion
   if tens_strength == 0 or np.isnan(tens_strength):
       tens_strength = cohesion/2
   if cr == 0 or np.isnan(cr):
       cr = def_cr/1e6
   if cw == 0 or np.isnan(cw):
       cw = def_cw/1e6
   if salinity == 0 or np.isnan(salinity):
       salinity = def_salinity/1e6
   if dens_c == 0 or np.isnan(dens_c):
       dens_c = np.float((eos(T0_mean, pres0_mean, salinity,0))[1])       
   if visc_c == 0 or np.isnan(visc_c):
       visc_c = (eos(T0_mean, pres0_mean, salinity, dens_c))[2]
   if visc_w == 0 or np.isnan(visc_w):
       visc_w = (eos(T0_mean, pres0_mean,salinity,0))[0]


#calculate some useful parameters
   s1 = s1_tot - pres0                                                          # effective maximum principal stress [MPa]
   s3 = stress_ratio*s1                                                         # effective minimum principal stress [MPa]
   theta = (1 - np.sin(np.radians(friction)))/(1 + np.sin(np.radians(friction)))
   p_lim_shear = (stress_ratio - theta)/(1 - theta)*s1 + cohesion*np.cos(np.radians(friction))/np.sin(np.radians(friction))
   p_lim_tensile = s3 + tens_strength
   p_lim = min(p_lim_shear, p_lim_tensile)                                      # limit overpressure, either for shear failure or tensile failure [MPa]

   gamma = visc_c/visc_w                                                        # non-dimensional viscosity ratio [-]
   delta = (visc_w-visc_c)/visc_w
   omega = (visc_c+visc_w)/(visc_c-visc_w)*np.log(np.sqrt(visc_c/visc_w))-1
   compr = cr+por*cw

   return thick, area_res, perm, por, dens_c, visc_c, visc_w, compr, p_lim, rc, gamma, delta, omega