# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 11:47:40 2016

@author: mbialoncz
"""

from __future__ import division
import numpy as np
import scipy.linalg as lin
import scipy.optimize as opt
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys
import random


from mean_field import *

mu_values = np.arange(0.4,0.7,0.001)

kappa = 1.
mu = 1.2
Q=0.2

Q_values = np.arange(-0.6, 0.6,0.005)

E_values = [Energy([q,0],[0,0],0,0,mu,1.,6,'K01') for q in Q_values]
plt.plot(Q_values, E_values)
plt.show()

n_values = [ParticleNumber([Q,0], [0,0] , 0, 0, mu, kappa, 6, 'K01') for mu in mu_values]
plt.plot(mu_values, n_values)
plt.show()

print Bis(0, 5, [Q,0], [0, 0], 0, 0, kappa, 6, 'K01')

E_values = []  
n_values = []
Q0 = [0.1]

for mu in mu_values : 
    
    f = lambda x : Energy([x[0],0],[0., 0.], 0, 0, mu, 1., 6, 'K01')
    minim = opt.minimize(f, Q0, method = 'Nelder-Mead', tol = 1e-6)
    print mu, minim.x[0], ParticleNumber([minim.x[0],0], [0,0], 0, 0, mu, kappa, 6, 'K01')
    n_values.append(ParticleNumber([minim.x[0],0], [0,0], 0, 0, mu, kappa, 6, 'K01'))
    E_values.append(minim.fun)
    Q0 = [minim.x[0]+0.1*mu/50]
    

plt.figure(1)
plt.subplot(211)
plt.plot(mu_values, E_values)

plt.subplot(212)
plt.plot(mu_values, n_values)
plt.show()
