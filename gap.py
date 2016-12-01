

import numpy as np
import scipy.linalg as lin
import scipy.optimize as opt
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys
import random


from mean_field import *

kappa=1.

#Q_values = np.arange(-3. ,3., 0.1)
#E_values = [Energia(0 , 5, [q,0],[0,0],0,0,kappa, 6, 'T1') for q in Q_values]

f = lambda x : Energia(0, 5, [x[0],0], [0,0],0,0,kappa, 6, 'T1')
Q0 = [0.2]
minim = opt.minimize(f, Q0, method='Nelder-Mead', tol=1e-6)

print minim.x
mu  = Bis(0,5,[minim.x[0],0],[0,0],0, 0, kappa, 6, 'T1' )


print 'spectral gap=' , spectral_gap([minim.x[0],0],[0,0],0,0, kappa, mu, 6, 'T1')

plt.plot(Q_values, E_values)
plt.show()


