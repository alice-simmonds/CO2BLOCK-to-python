# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 09:36:53 2022

@author: simmonds
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


fname = pd.read_excel('example_data.xlsx')   #name of input data file, no directory required if saved in the same location as the python functions
data = pd.DataFrame(fname)

data.head

correction = 'off'
dist_min = 0.5
dist_max = 'auto'
nr_dist = 15
nr_well_max = 20
rw = 0.156
time_yr = 30
maxQ = 2


from calculate import calculate

d_list,well_list,d_max,Q_M_each,V_M,Table_Q,Table_V = calculate(data,correction,dist_min,dist_max,nr_dist,nr_well_max,rw,time_yr,maxQ)



# sustainable flow rate per well
fig, ax1 = plt.subplots(figsize=(4,3),dpi=100)
Qcontour = ax1.contour(d_list,well_list,np.real(Q_M_each),10,colors='grey',linewidths=1)
ax1.clabel(Qcontour,fmt="%.1f")
ax1.set_title('Maximum sustainable per well flow rate $Q_{M}$ (Mt/yr)',fontsize=15)

plt.plot(d_max,well_list,linewidth=2, color='red')
ax1.set_xlim([min(d_list),max(d_list)])
ax1.set_xlabel('inter-well distance $d$ (km)',fontsize=10)
ax1.set_ylabel('number of wells $n$',fontsize=10)



plt.show()


# sustainable storage capacity
fig, ax2 = plt.subplots(figsize=(4,3),dpi=100)
Vcontour = ax2.contour(d_list,well_list,np.real(V_M),24,colors='grey',linewidths=1)
ax2.clabel(Vcontour,fmt="%.0f")
ax2.set_title('Maximum sustainable storage $V_{M}$ (Mt)',fontsize=15)

plt.plot(d_max,well_list,linewidth=2, color='red')
ax2.set_xlim([min(d_list),max(d_list)])
ax2.set_xlabel('inter-well distance $d$ (km)',fontsize=10)
ax2.set_ylabel('number of wells $n$',fontsize=10)



plt.show()
