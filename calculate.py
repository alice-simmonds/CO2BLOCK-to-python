# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 15:39:44 2022

@author: simmonds
"""


import pandas as pd
import numpy as np
from scipy.special import lambertw


def calculate(data,correction,dist_min,dist_max,nr_dist,nr_well_max,rw,time_yr,maxQ):
    
    #read data
    from read_data import read_data
    thick, area_res, perm, por, dens_c, visc_c, visc_w, compr, p_lim, rc, gamma, delta, omega =  read_data(data)
    time = time_yr*86400*365                                                # injection time [seconds]
    R_influence = np.sqrt(2.246*perm*time/(visc_w*compr))                   # pressure propagation radius for the time of injection

    
    if nr_well_max == 'auto':                                               # calculate maximum well number if not set
        nr_well_max = np.floor(area_res/(dist_min**2))
      
    if dist_max == 'auto':                                                  # calculate maximum interwell distance if not set
        dist_max = np.sqrt(2*area_res)/2
    
    d_list = np.linspace(dist_min, dist_max, nr_dist)                       # interwel distance list
    M0 = perm/1e-13                                                         # guess value for total injection rate [Mton/y]
    

    well_list = []
    d_max = []
    w_id = 0
    
    Q0_vec = []
    d_min_p = []
    w_id_max = 0
    
    
    # loop to give dimensions of array for superposed pressure
    for x_grid_num2 in range(1,np.int_(np.sqrt(nr_well_max))+1):           
        if x_grid_num2*(x_grid_num2+1) < nr_well_max:
            plus2 = 1
        else:
            plus2 = 0
        
        for y_grid_num2 in  range(x_grid_num2,x_grid_num2+plus2+1):
            w_id_max = w_id_max + 1
            p_sup_vec = np.ones([w_id_max,nr_dist])
            
    
    from Nordbotten_solution import Nordbotten_solution

    
    for x_grid_num in range(1,np.int_(np.sqrt(nr_well_max))+1):             # number of wells on a horizontal row
        if x_grid_num*(x_grid_num+1) < nr_well_max:
            plus = 1
        else:
            plus = 0
        
        
        for y_grid_num in range(x_grid_num,x_grid_num+plus+1):              # number of wells on a vertical row
            w_id = w_id + 1                                                 # well number scenario ID
            w = x_grid_num*y_grid_num                                       # well number for each scenario
            well_list.insert(w_id,w)                                        # store in list
            d_max.insert(w_id, np.sqrt(area_res/w))                         # maximum interwell distance for each well number [km]
           
            #calculate pressure build-up for a reference flow rate Q0 at each scenario, moved these out of the nested for loop
            Q0 = M0*1e9/dens_c/365/86400/w                                  # injection rate per well [m3/s]
            Q0_vec.insert(w_id,Q0)                                          # store in list
            csi = np.sqrt(Q0*time/np.pi/por/thick)                          # average plume extension [m]
            psi = np.exp(omega)*csi                                         # equivalent plume extension [m]
            p_c = (Q0*visc_w)/(2*np.pi*thick*perm)/1e6                      # characteristic pressure [MPa]
            d_min_p.insert(w_id,2*csi/1000)                                 # minimum interwell distance for each well number case [km]

            

            for d in range(1,nr_dist+1):
                #calculate grid and distances
                distance = d_list[d-1]*1000
             
                wells_coord_x = np.tile(np.arange(0,distance*x_grid_num-1,step=distance),(y_grid_num,1))
                wells_coord_y = np.tile(np.matrix.transpose(np.array([np.arange(0,distance*y_grid_num-1,step=distance)])),(1,x_grid_num))
                central_well_x = np.int_(np.ceil(x_grid_num/2))                                     # position of central well in x coord
                central_well_y = np.int_(np.ceil(y_grid_num/2))                                     # position of central well in y coord
                        
                dist_vec_x = wells_coord_x - wells_coord_x[central_well_y-1,central_well_x-1]       # distance in x from central well [km]
                dist_vec_y = wells_coord_y - wells_coord_y[central_well_y-1,central_well_x-1]       # distance in y from central well [km]
            
                dist_vec = np.sqrt(np.square(dist_vec_x) + np.square(dist_vec_y))                   # distance from central well
                dist_vec[central_well_y-1,central_well_x-1] = rw                                    # assign well's radius to the central well
    

                p_sup = 0
                for i in range(1,w+1):
                    dist_vec_flat = np.ndarray.flatten(dist_vec,order='F')
                    r = np.float(dist_vec_flat[i-1])
                    Delta_p = Nordbotten_solution(r,R_influence,psi,rc,gamma)*p_c                   # overpressure according to Nordbotten and Celia solution for overpressure [MPa]
                    p_sup = p_sup + Delta_p                                                         # superposed overpressure [MPa]
                    
                b = np.ones([1,w_id])
                if correction == 'off':
                    sup_error = 0
                  
                    b[:] = (visc_w-visc_c)/4/np.pi/perm/thick
                elif correction == 'on':                                                            # correction for superposition error (De Simone et al.)
                    if w < 9 or R_influence*csi/(distance**2) < 1:
                        sup_error = 0
                        
                        b[:] = (visc_w-visc_c)/4/np.pi/perm/thick
                    else:
                        sup_error = w*delta/4 * np.log(R_influence*csi/(distance**2))
                        b[-1] = ((visc_w-visc_c)/4/np.pi/perm/thick*(1+w/4))
                p_sup = p_sup - sup_error*p_c
 
                p_sup_vec[w_id-1,d-1] = p_sup
    
 
    #calculate injectable per well flow rate Q_M_each and a total storage capacity V_M for each scenario
    b = np.tile(np.reshape(b,(w_id,1)), (1, nr_dist))
    q1 = np.tile(np.reshape(Q0_vec,(w_id,1)),(1,nr_dist))
    well_mat = np.tile(np.reshape(well_list,(w_id,1)), (1, nr_dist))
    

    p1 = p_sup_vec*np.int(1e6)
    p2 = p_lim*1e6

  
    q2 = - p2/b/(lambertw(-p2/q1/b*np.exp(-p1/q1/b),-1))                                            # [m3/s] limit flow rate at each well with non-linear multi-phase relationship
    Q_M_each = q2*86400*365*dens_c/1e9                                                              # [Mton/year]
  

    
    #rescale according to lower constraints
    for dd in range(1,nr_dist+1):
        distance = d_list[dd-1]*1000
        if Q_M_each[0,dd-1] > 0.9999*maxQ:              # rescale for n_well = 1
            Q_M_each[0,dd-1] = 0.9999*maxQ
        for nn in range(2,w_id+1):                      # rescale for n_well > 1
            if Q_M_each[nn-1,dd-1] > 0.9999*maxQ or q2[nn-1,dd-1] > distance**2*np.pi*por*thick/4.001/time:
                Q_M_each[nn-1,dd-1] = min(0.9999*maxQ , distance**2*np.pi*por*thick/4.001/time*86400*365*dens_c/1e9)
    
    Q_M_tot = Q_M_each*well_mat                         # [Mton/year] total sustainable flow rate
    V_M = Q_M_tot*time_yr                          # [Mton] total sustainable injectable mass
    
    
    #upper constraint
    d_mat = np.tile(d_list,[w_id,1])
    d_max = np.tile(d_max,[nr_dist,1])
    d_max = np.matrix.transpose(d_max)
    d_max_check = d_mat/d_max
    

    #apply upper constraint
    possible = d_max_check < 1
    Q_poss_each = np.real(Q_M_each * possible)
    V_poss = np.real(V_M * possible)
    
    #write results in tables
    varNamesQ = ['Q_M_for_d_'+str(np.int(d_list[s-1]*1000))+'_m' for s in range(1,nr_dist+1)]
    Table_Q = pd.DataFrame(Q_poss_each, columns= varNamesQ, index= well_list)
            
    varNamesV = ['V_M_for_d_'+str(np.int(d_list[s-1]*1000))+'_m' for s in range(1,nr_dist+1)]
    Table_V = pd.DataFrame(V_poss, columns= varNamesV, index= well_list)  
  
    Table_Q.to_excel("Q_M_max_per_well_inj_rate.xlsx",index_label = "number_of_wells")    
    Table_V.to_excel("V_M_max_storage_capacity.xlsx",index_label = "number_of_wells")   
    
    
    return d_list,well_list,d_max,Q_M_each,V_M,Table_Q,Table_V