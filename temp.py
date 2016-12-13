# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
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

Q_values = np.arange(-3,3,0.5)
E_values = [Energia(0,7, [q,0], [0,0], 0, 0,  1., 6, 'K01') for q in Q_values] 
plt.plot(Q_values, E_values)
plt.show()