# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 10:55:11 2016

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

#file contains lines with the data : kappa, Q_minimal, mu_minimal, E_minimal
with open('result_K01_6', 'r') as res :
    results = []
    for line in res : 
        print line
        results.append([float(x) for x in line.split()])
        
kappa = results[0][0]
Q_minimal = results[0][1]
mu_minimal = results[0][2]
print dispersion(0.3, 0.5, [Q_minimal,0], [0,0], 0, 0, mu_minimal, kappa, 6, 'K01')

#computation of spectral gap
#print spectral_gap([minim.x[0],0],[0,0],0,0,kappa, mu, 6, 'T1')

#partial plots of the brillouin zone
#kx, ky = 0.5,2.0933
#k_values = np.linspace(-2*math.pi, 2*math.pi, 1000)
#omega_1 = lambda k1 : dispersion(k1,ky, [Q_minimal, 0],[0,0],0, 0, mu_minimal, kappa, 6, 'T1')[0]
#omega_1_values = [omega_1(x) for x in k_values]
#plt.plot(k_values, omega_1_values)
#plt.show()
print modes_of_dispersion([Q_minimal,0], [0,0], 0, 0, mu_minimal, kappa, 6, 'K01')