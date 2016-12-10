# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 12:51:25 2016

@author: mbialoncz
"""

import numpy as np
import scipy.linalg as lin
import scipy.optimize as opt
import math
import random
import matplotlib.pyplot as plt
import sys
import random as rand

from mpl_toolkits.mplot3d import Axes3D
from mean_field import *
from pylab import *

Ns = 12
with open('result_T1_12', 'r') as res :
    results = []
    for line in res : 
        print line
        results.append([float(x) for x in line.split()])
        
kappa       = results[5][0]
Q_minimal   = results[5][1]
mu_minimal  = results[5][2]

H_values = np.arange(0, 2, 0.01)
magnetization_values = []

for H in H_values :
    #find the wave vector corresponding to zero mode of energy
    g = lambda k : np.abs(dispersion(k[0],
                                     k[1],
                                     [Q_minimal,0],
                                     [0,0],
                                     0,
                                     0,
                                     mu_minimal,
                                     kappa,
                                     Ns,
                                     'T1')[0] - H/2)
    sol = opt.minimize(g, [-2.,2.], method = 'Nelder-Mead')
    k_minim = sol.x
    n1 = Ns * k_minim[0]/(2 * math.pi)
    n2 = Ns * k_minim[1]/(2 * math.pi)

    M, dim = HamiltonianMatrix(n1,
                               n2,
                               [Q_minimal,0],
                               [0,0],
                               0,
                               0,
                               mu_minimal,
                               kappa,
                               Ns,
                               'T1')
    B = np.identity(dim)
    B[dim/2:dim, dim/2:dim] = -np.identity(dim/2)

    w, v = lin.eig(np.dot(B, M))
    # print 'condensate:', v[0]

    # normalization factor
    r = np.abs(v[0][0])**2 + np.abs(v[0][1])**2      

    result = (np.abs(v[0][0])**2 - np.abs(v[0][1])**2)/r
    magnetization_values.append(result)

    #only one condensate - computing z component of magnetization?
    print 'Result:', result 
    
plt.plot(H_values, magnetization_values)
plt.xlabel('H')
plt.ylabel('S_z')
plt.savefig('magnetization_triangular')
plt.show()
