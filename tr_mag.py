# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 11:31:20 2016

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



kappa = 1.

H_values = np.arange(0,1,0.1)
M_values = []

for H in np.arange(0.0,3.1,0.05) :
        
        f = lambda x: Energia(0,5, x[0], 0, 0, H, kappa,6,"Sa1")
        x00 = 0.2;
        x0 = np.zeros(1,dtype = float)
        x0[0] = x00
        minim = opt.minimize(f, x0, method='Nelder-Mead')
        Q0 = (minim.x)[0]
        
        mu = Bis(0., 5., Q0, 0, 0, H, kappa, 6, 'Sa1')
        
#        M = AverageMagnetization(Q0, 0,0, H, mu, kappa, 6, 'Sa1')
#        
#        with open('MagnetizationSa1','a') as f :
#              f.write(str(H) + " " + str(M)+'\n')
#        
#        print "wartosci", H, M
        
        G = spectral_gap(Q0,0,0, H,  kappa, mu, 6, 'Sa1')
        with open('SpectralGapSa1','a') as f :
              f.write(str(H) + " " + str(G)+'\n')
       
        M_values.append(M)

#plt.xlabel('H')
#plt.ylabel('M')
#plt.plot(H_values, np.array(M_values))
#plt.savefig('MagnetizationSa1')

  
        
        