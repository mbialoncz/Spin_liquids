# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 12:55:25 2016

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

#kappa = 1.
#H=100
#
#f = lambda x: Energia(0,5, x, 0, 0, H, kappa,6,"Sa1") 
#x0 = np.zeros(2, dtype = complex)
#
#x0[0] = 0.2
#x0[1] = 0.5
#
#minim = opt.minimize(f, x0, method='Nelder-Mead')
#
#print minim.x
#print minim.f
BogolubovTransformation(2,5 , [0.2,1], [0.45,1], 1., 0, 1.5, 1., 6, 'K01') 