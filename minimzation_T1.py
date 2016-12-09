import sys
import math
import random
import numpy as np
import random as rand
import scipy.linalg as lin
import scipy.optimize as opt
import matplotlib.pyplot as plt

from mean_field import *
from pylab import *

def find_q_mu(kappa, Ns):
    """ Such that energy is minimized """
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

    return Q_minimal, mu_minimal

def minimize_T1():
    """ Find Q and mu paramters for multiple kappas """
    # Set some paremeters
    Ns = 12
    # Prepare empty file
    savepath = 'results/T1_{}.dat'.format(Ns)
    open(savepath, 'w').close()

    # Declare research space
    kappas = np.arange(1.2, 1.6, 0.2)
    for kappa in kappas:
        Q, mu = find_q_mu(kappa, Ns)

        # Save results at every iteration
        scoreline = '{} {} {}\n'.format(kappa, Q, mu)
        with open(savepath, 'a') as fout:
            fout.write(scoreline)

    return savepath

