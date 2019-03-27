#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 15:27:32 2019

@author: btreece
"""

import numpy as np
from matplotlib import pyplot as plt

workingdir = "/home/btreece/Programs/BASIC_MPI/"

filer = open(workingdir + "OUTPUT.txt")
lines = filer.readlines()
filer.close()

#x  = []
#v  = []
#xx = []
#vv = []
#xv = []

for i in range(len(lines)):
    if lines[i][:7] == "NEW_RUN":
        i += 1
        tmp = lines[i].split()
        dt = float(tmp[0])
        dT_out = float(tmp[1])
        m = float(tmp[2])
        eta = float(tmp[3])
        kT = float(tmp[4])
        F = float(tmp[5])
        n_samples = int(tmp[6])
        while (lines[i][:7] != "NEW_RUN"  and i < len(lines)-1):
            i+=1
            if lines[i][:4] == "<x>:":
                tmp = lines[i].split()[1:]
                x = [float(j) for j in tmp]
            if lines[i][:4] == "<v>:":
                tmp = lines[i].split()[1:]
                v = [float(j) for j in tmp]
            if lines[i][:5] == "<xx>:":
                tmp = lines[i].split()[1:]
                xx = [float(j) for j in tmp]
            if lines[i][:5] == "<vv>:":
                tmp = lines[i].split()[1:]
                vv = [float(j) for j in tmp]
            if lines[i][:5] == "<xv>:":
                tmp = lines[i].split()[1:]
                xv = [float(j) for j in tmp]
                
def brownian_msd(t, kT, eta, m):
    return (2*kT/eta)*(t - (2*m/eta)*(1 - np.exp(-eta*t/m)) + (0.5*m/eta)*(1 - np.exp(-2.0*eta*t/m)))

def brownian_msv(t, kT, eta, m):
    return (kT/m)*(1 - np.exp(-2.0*eta*t/m))

def brownian_cov_x_v(t, kT, eta, m):
    return (kT/eta)*(1-np.exp(-eta*t/m))**2.0