import numpy as np
import scipy.linalg as lin
import scipy.optimize as opt
import math
import random
import matplotlib.pyplot as plt
import sys
import random as rand

from mean_field import *
from pylab import *

# generates file 'result_T11_12' with the following structure:
# kappa  Q_minimal   mu_minimal  E_minimal
# where Q_minimal , mu_minimal are parameters
# for which there is minimum of energy for a given kappa an

#plot
#Q_values = np.arange(-3. ,3., 0.1)
#E_values = [Energia(0 , 5, [q,0],[0,0],0,0,kappa, 6, 'T1') for q in Q_values]

Ns = 12

with open('result_T1_12', 'wr') as res: 
    for kappa in np.arange(1.2, 1.6, 0.2) :
        print 'Minimizing energy for kappa:', kappa
        f = lambda x : Energia(0,
                               5,
                               [x[0], 0],
                               [0, 0],
                               0,
                               0,
                               kappa,
                               Ns,
                               'T1')
        Q0 = [0.2]
        minim = opt.minimize(f, Q0, method='Nelder-Mead', tol=1e-6)
        Q_minimal = minim.x[0]

        print 'Found Q:', Q_minimal

        mu_minimal  = Bis(0,
                          5,
                          [Q_minimal, 0],
                          [0, 0],
                          0,
                          0,
                          kappa,
                          Ns,
                          'T1')

        print 'Found mu:', mu_minimal

        # Save to file
        scoreline = '{} {} {} {}\n'.format(kappa,
                                           Q_minimal,
                                           mu_minimal,
                                           minim.fun)
        res.write(scoreline)
