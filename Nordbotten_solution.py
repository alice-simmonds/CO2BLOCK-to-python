# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 16:45:19 2022

@author: simmonds
"""
import numpy as np

from FD_Nor import FD_Nor

def Nordbotten_solution(r,R_influence,psi,rc,gamma):
    if r <= psi:
        PD = gamma*np.log(psi/r) + FD_Nor(psi,R_influence,rc)
    elif r > psi and r <= R_influence:
        PD = FD_Nor(r,R_influence,rc)
    else:
        PD = 0

    return PD