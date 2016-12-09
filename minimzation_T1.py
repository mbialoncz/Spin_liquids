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

#generates file 'result_T11_12' with the following structure:
#kappa  Q_minimal   mu_minimal  E_minimal
#where Q_minimal , mu_minimal are parameters for which there is minimum of energy for a given kappa an


#plot
#Q_values = np.arange(-3. ,3., 0.1)
#E_values = [Energia(0 , 5, [q,0],[0,0],0,0,kappa, 6, 'T1') for q in Q_values]


Ns = 12


with open('result_T1_12', 'wr') as res : 
    for kappa in np.arange(1.2, 1.6, 0.2) :
        f = lambda x : Energia(0, 5, [x[0],0], [0,0],0,0,kappa, Ns, 'T1')
        Q0 = [0.2]
        minim = opt.minimize(f, Q0, method='Nelder-Mead', tol=1e-6)
        Q_minimal = minim.x[0]
        mu_minimal  = Bis(0,5,[Q_minimal,0],[0,0],0, 0, kappa, Ns, 'T1')
    
        print kappa, Q_minimal, minim.fun
        res.write(str(kappa) + ' ' + str(Q_minimal) +' ' +str(mu_minimal)+ ' '+ str(minim.fun)+'\n')




#computation of spectral gap
#print spectral_gap([minim.x[0],0],[0,0],0,0,kappa, mu, 6, 'T1')

#partial plots of the brillouin zone



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


