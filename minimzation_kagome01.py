# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 15:48:16 2016

@author: mbialoncz
"""
#generates file 'result_K01_6' with the following structure:
#kappa  Q_minimal   mu_minimal  E_minimal
#where Q_minimal , mu_minimal are parameters for which there is minimum of energy for a given kappa and ansatz

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


with open('result_K01_12_v2', 'wr') as res : 
    for kappa in np.arange(0.5, 2., 0.5) :
        f = lambda x : Energia(0, 5, [x[0],0], [0,0],0,0,kappa, Ns, 'K01')
        Q0 = [0.2]
        minim = opt.minimize(f, Q0, method='Nelder-Mead', tol=1e-6)
        Q_minimal = minim.x[0]
        mu_minimal  = Bis(0,5,[Q_minimal,0],[0,0],0, 0, kappa, Ns, 'K01')
    
        print kappa, Q_minimal, minim.fun
        res.write(str(kappa) + ' ' + str(Q_minimal) +' ' +str(mu_minimal)+ ' '+ str(minim.fun)+'\n')