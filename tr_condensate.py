# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 11:37:13 2017

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

N = 6
kappa = 1.

def Energy_condensate_full(Q, F1, x, y, H, mu, kappa, Ns) : 
    
    result = 0
    for n1 in range(Ns) : 
        for n2 in range(Ns) : 
            M, dim = HamiltonianMatrix(n1, n2, Q, F1, 0, H, mu, kappa, Ns, 'T1') 
        
            B = np.identity(dim)
            B[dim/2:dim, dim/2:dim] = -np.identity(dim/2)
        
            eig = np.absolute(np.real(lin.eigvals(np.dot(B,M))))
        
            result += sum(eig)/2
            
            vec = [x[Ns * n1 + n2], np.conjugate(y[Ns * ((Ns - n1) % Ns) + (Ns - n2) % Ns])]
            
            result += np.dot(vec, np.dot(np.conj(vec).T, M))
            
    return result + 3 * Ns ** 2 * (np.abs(F1)**2 - np.abs(Q)**2) - Ns * (mu + kappa) + Ns * H
    

 
def Average_number(Q, F1, x, y, H, mu, kappa, Ns) :
     
     result = 0
     
     for n1 in range(Ns) : 
        for n2 in range(Ns) : 
           result += np.abs(x[Ns * n1 + n2])**2 + np.abs(y[Ns * n1 + n2])**2
           
     return ParticleNumber(Q, F1, 0, H, mu, kappa, Ns, 'T1') + result/(Ns*Ns)



def Energy_condesate(params, H, Ns, kappa) :
    Q = params[0]
    F1 = params[1]
    mu = params[2]
    xc1 = params[3] + params[4] * 1j
    yc1 = params[5] + params[6] * 1j
    xc2 = params[7] + params[8] * 1j
    yc2 = params[9] + params[10] * 1j   
    
    x = np.zeros(Ns**2)
    y - np.zeros(Ns**2)
    x[Ns**2/3 + Ns/3] = xc1
    y[Ns**2/3 + Ns/3] = yc1
    x[2/3 * Ns**2 + 2/3 * Ns] = xc2
    y[2/3 * Ns**2 + 2/3 * Ns] = yc2
    
        
    
    return Energy_condensate_general1(Q, F1, x, y, H, mu, kappa, Ns)
    
def Bis_condensate(a, b, Q, F1, x, y, H, kappa, Ns) :
    
#     while (ParticleNumber(Q, F1, F2, H, a, kappa, Ns, ansatz) == np.inf or 
#             ParticleNumber(Q, F1, F2, H, a, kappa, Ns, ansatz) == np.nan)  : 
#        
#        while ParticleNumber(Q, F1, F2, H, b, kappa, Ns, ansatz) < kappa :
#            b=(a+b)/2
#        
#        a = (a+b)/2
#    
#    if ParticleNumber(Q, F1, F2, H, a, kappa, Ns, ansatz) <= kappa : 
#        a = 2*a - b
    c = 0 
    it = 0

    while (Average_number(Q, F1, x, y,  H, b, kappa,Ns) == np.inf or
            np.isnan(Average_number(Q, F1, x,y, H, b, kappa,Ns)) or
              Average_number(Q, F1,x, y, H, b, kappa,Ns) >= kappa):
                  b = 2 * b
    
    while (Average_number(Q, F1,x,y, H, a, kappa,Ns) == np.inf or 
            np.isnan(Average_number(Q, F1,x, y, H, a, kappa,Ns))  or 
                np.absolute(Average_number(Q, F1, x, y, H, a, kappa,Ns) -  
                      Average_number(Q,F1, x, y, H, b, kappa, Ns)) > 0.00001) :
        
        c = (a+b)/2
        it += 1
        
        if it > 30:
            return c
        
        #print Q, c
        
        if (np.sign(Average_number(Q, F1, x, y, H, c, kappa, Ns) - kappa) == 
                np.sign(Average_number(Q, F1, x, y, H, b, kappa, Ns)-kappa)) :
            b = c
        else : 
            a = c
            
    return c

def Energy_condensate_reduced(Q, F1, x, y, H, kappa, Ns) :
    mu = Bis_condensate(0, 10, Q, F1, x, y, H, kappa, Ns)
    return Energy_condensate_full(Q, F1, x, y, H,mu,  kappa, Ns)
    

x0 = np.zeros(36, dtype = complex)
y0 = np.zeros(36, dtype = complex)

for n1 in range(N) :
    for n2 in range(N) :
        x0[N * n1 + n2] = 0
        y0[N * n1 + n2] = 0

#x0[N**2/3 + N/3] = np.random.rand() + np.random.rand()*1j
#y0[N**2/3 + N/3] = np.random.rand() + np.random.rand()*1j

Q0 = 0.2
F10 = 1.
mu0 = 1.2   
x0[N**2/3 + N/3] = 1. 
y0[N**2/3 + N/3] = 1j
x0[2/3 * N**2 + 2/3 * N] = 1
y0[2/3 * N**2 + 2/3 * N] = -1j 

f = lambda mu : Energy_condensate_general(Q0, F10, x0, y0, H, mu , kappa, N) 
 
a = Bis_condensate(0, 5, Q0, F10, x0, y0, 0. , kappa, N)
print Average_number(Q0, F10, x0, y0, 0., a, kappa, N) 

mu_values = np.arange(5.,10.,0.1)

Q_values = np.arange(-5,5, 1.)
F_values = np.arange(-0.5,0.5, 0.1)
particle_values =  [Average_number(Q0, F10, x0, y0, 0., mu, kappa, N) for mu in mu_values]
Energy_values = [np.abs(Energy_condensate_reduced(Q0, F, x0, y0, 0.0, kappa, N)) for F in F_values]
#print particle_values
print Energy_values



plt.plot(Q_values, Energy_values)
plt.show() 
     

     
     
           
    
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
 
            
    