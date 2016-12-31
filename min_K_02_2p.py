# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 10:54:04 2016

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


with open('min_K_02_2p', 'wr') as res : 
    for kappa in np.arange(0.5, 2., 0.5) :
        f = lambda x : Energia(0, 5, [x[0],0], [x[1],0],0,0,kappa, Ns, 'K02')
        x0 = [0.2,0.2]
        minim = opt.minimize(f, x0, method='Nelder-Mead', tol=1e-6)
        Q_minimal = minim.x[0]
        F1_minimal = minim.x[1]
        mu_minimal  = Bis(0,5,[Q_minimal,0],[F1_minimal,0],0, 0, kappa, Ns, 'K02')
    
        print kappa, Q_minimal, F2_minimal, minim.fun
        res.write(str(kappa) + ' ' + str(Q_minimal) +' ' +str(F1_minimal)+' '+str(mu_minimal)+ ' '+ str(minim.fun)+'\n')