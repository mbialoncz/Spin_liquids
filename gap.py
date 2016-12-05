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

kappa=1.25

#plot
#Q_values = np.arange(-3. ,3., 0.1)
#E_values = [Energia(0 , 5, [q,0],[0,0],0,0,kappa, 6, 'T1') for q in Q_values]

#minimization
f = lambda x : Energia(0, 5, [x[0],0], [0,0],0,0,kappa, 6, 'T1')
Q0 = [0.2]
minim = opt.minimize(f, Q0, method='Nelder-Mead', tol=1e-6)
Q_minimal = minim.x[0]
mu  = Bis(0,5,[Q_minimal,0],[0,0],0, 0, kappa, 6, 'T1')

with open('result', 'wr') as res : 
    res.write(kappa, ' ', Q_minimal, ' ', minim.fun)

#computation of spectral gap

print Energia(0,5, [minim.x[0],0], [0,0],0,0,kappa, 6 ,'T1')
print minim.x[0]
#print spectral_gap([minim.x[0],0],[0,0],0,0,kappa, mu, 6, 'T1')

#partial plots of the brilluin zone
kx, ky = 0.5,0.5
k_values = np.linspace(-2*math.pi, 2*math.pi, 100)
omega_1 = lambda k1 : dispersion(k1,ky, [Q_minimal, 0],[0,0],0, 0, mu, kappa, 6, 'T1')
omega_1_values = [omega_1(x) for x in k_values]
plt.plot(k_values, omega_1_values)
plt.show()



#g = lambda k : dispersion_0(k[0], k[1], [Q_minimal, 0],[0,0],0, 0, mu, kappa, 6, 'T1')
#
#list_of_minimas = []
#
#for _ in xrange(3) :
#   k0 = [2*math.pi * rand.random(), 2 * math.pi * rand.random()]    
#   k_minim = opt.minimize(g, k0, method = 'Nelder-Mead', tol = 1e-6)
#   list_of_minimas.append([k_minim.x,k_minim.func])
#   
#print np.sort(list_of_minimas)[0]
#
#
#
#   
#
#with open('results', 'wr') as r :
#    r.write(')
#
#print k_minim.x
#print k_minim.fun




#plt.plot(Q_values, E_values)
#plt.show()


