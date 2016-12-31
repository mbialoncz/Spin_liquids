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
with open('result_K01_12_v2', 'r') as res :
    results = []
    for line in res : 
        print line
        results.append([float(x) for x in line.split()])
        
kappa = results[1][0]
Q_minimal = results[1][1]
mu_minimal = results[1][2]
print 'kappa=', kappa

#print modes_of_dispersion([Q_minimal,0], [0,0], 0, 0, mu_minimal, kappa, 12, 'T1')
#print dispersion(-2*math.pi/10, -2*math.pi/10,[Q_minimal,0], [0,0], 0, 0, mu_minimal, kappa, 6, 'T1')
#print modes_of_dispersion([Q_minimal,0], [0,0], 0, 0, mu_minimal, kappa, 6, 'T1')

#print modes_of_dispersion([Q_minimal,0], [0,0], 0, 0, mu_minimal, kappa, Ns, 'T1')

H_values = np.arange(0,10,0.1)
magnetization_values = []

for H in H_values :
    #find the wave vector corresponding to zero mode of energy
    g = lambda k : np.abs(dispersion(k[0], k[1], [Q_minimal,0], [0,0], 0, 0, mu_minimal, kappa, Ns, 'K01')[0] - H/2)
    MINIMAL = 10000
 #   for _ in range(10) :
 #   k1 = random.random()
    sol = opt.minimize(g, [-2.,2.], method = 'Nelder-Mead')
    k_minim = sol.x
    print "solucja", k_minim/math.pi, sol.fun, g([0,0])
#    print sol.x
#    print sol.fun
    n1 = Ns * k_minim[0]/(2 * math.pi)
    n2 = Ns * k_minim[1]/(2 * math.pi)
    M, dim = HamiltonianMatrix(n1, n2, [Q_minimal,0], [0,0], 0, 0, mu_minimal, kappa, Ns, 'K01')
    B = np.identity(dim)
    B[dim/2:dim, dim/2:dim] = -np.identity(dim/2)
    w, v = lin.eig(np.dot(B,M))
    print "eigenvalues", w-H/2
    c = v[2]
    print c                                        #condensate
    r = np.sum(np.abs(c)**2)    #normalization
    m = np.sum(np.abs(c[:3])**2) - np.sum(np.abs(c[3:])**2)
    magnetization_values.append(m/r)
    #only one condensate - computing z component of magnetization?
    print 'magnetization', m
    
plt.plot(H_values, magnetization_values)
plt.xlabel('H')
plt.ylabel('S_z')
plt.savefig('magnetization_kagomer')
plt.show()

print modes_of_dispersion([Q_minimal,0], [0,0], 0, 0, mu_minimal, kappa, Ns, 'K01')
    
    


#kx, ky = 0.5,2.0933
#k_values = np.linspace(-2*math.pi, 2*math.pi, 1000)
#omega_1 = lambda k1 : dispersion(k1,ky, [Q_minimal, 0],[0,0],0, 0, mu_minimal, kappa, 6, 'T1')[0]
#omega_1_values = [omega_1(x) for x in k_values]
#plt.plot(k_values, omega_1_values)
#plt.show()
