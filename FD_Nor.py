# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 17:15:36 2022

@author: simmonds
"""

import numpy as np

def FD_Nor(x,R_influence,rc):
    if  x < rc:
        if R_influence <= rc:
            FD = np.log(R_influence/x)
        else:
            FD = np.log(rc/x) + 2/2.25*(R_influence/rc)**2 - 3/4
    else:
        FD = 0
        
    return FD