# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 22:15:00 2017

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
kappa = 1.3

def Energy_condensate_full(Q, F1, x, y, H, mu, kappa, Ns) : 
    
    if Q==0 and F1 ==0 :
        return 1e14
    
    m = find_minimum(Q, F1, mu, kappa, Ns)
    if m[0] < H/2 :
        return 1e14
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
            
    return result - 3 * Ns ** 2 * (np.abs(F1)**2 - np.abs(Q)**2)/2 - Ns**2 * mu*(1. + kappa) + Ns * H
    
def find_minimum(Q, F1, mu, kappa, Ns) :
    m = 100000000
    k_min = [0,0]
    for n1 in range(Ns):
        for n2 in range(Ns) : 
             k1 = 2 * n1 * np.pi/Ns
             k2 = 2 * n2 * np.pi/Ns
             a = mu + 2 * J * F1 * (np.cos(k1) + np.cos(k2) + np.cos(k1+k2))
             b = 2 * J * Q * (np.sin(k1) + np.sin(k2) - np.sin(k1+k2))
             Ek2 = a**2 - b**2
             if Ek2 < m :
                 m = Ek2
                 k_min = [k1,k2]
    return [m]+k_min
    

 
def Average_number(Q, F1, x, y, H, mu, kappa, Ns) :
     
     result = 0
     
     for n1 in range(Ns) : 
        for n2 in range(Ns) : 
           result += np.abs(x[Ns * n1 + n2])**2 + np.abs(y[Ns * n1 + n2])**2
           
     return ParticleNumber(Q, F1, 0, H, mu, kappa, Ns, 'T1') + result/(Ns*Ns)



def Energy_condensate(params, H, kappa, Ns) :
    Q = params[0]
    F1 = params[1]
#    F1 = params[1]
#    xc1 = params[2] + params[3] * 1j
#    yc1 = params[4] + params[5] * 1j
#    xc2 = params[6] + params[7] * 1j
#    yc2 = params[8] + params[9] * 1j   
    
    x = np.zeros(Ns**2)
    y = np.zeros(Ns**2)
#    x[Ns**2/3 + Ns/3] = xc1
#    y[Ns**2/3 + Ns/3] = yc1
#    x[2/3 * Ns**2 + 2/3 * Ns] = xc2
#    y[2/3 * Ns**2 + 2/3 * Ns] = yc2
    
        
    
    return Energy_condensate_reduced(Q, F1, x, y, H, kappa, Ns)
    
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
        
        if (np.sign(Average_number1(Q, F1, x, y, H, c, kappa, Ns) - kappa) == 
                np.sign(Average_number1(Q, F1, x, y, H, b, kappa, Ns)-kappa)) :
            b = c
        else : 
            a = c
            
    return c

def Energy_condensate_reduced(Q, F1, x, y, H, kappa, Ns) :

    mu = Bis_condensate(0., 5., Q, F1, x, y, H, kappa, Ns)

    return Energy_condensate_full(Q, F1, x, y, H,mu,  kappa, Ns)
    

def Average_number1(Q, F1, x, y, H, mu, kappa, Ns) :
     

     result = 0
     for n1 in range(Ns) :
         for n2 in range(Ns) :
             k1 = 2 * n1 * np.pi/Ns
             k2 = 2 * n2 * np.pi/Ns
             a = mu + 2 * J * F1 * (np.cos(k1) + np.cos(k2) + np.cos(k1+k2))
             b = 2 * J * Q * (np.sin(k1) + np.sin(k2) - np.sin(k1+k2))
#             print a/np.sqrt(a**2 - b**2)
             result += a/np.sqrt(a**2 -b**2) + np.abs(x[n1 * Ns + n2])**2 + np.abs(y[n1*Ns + n2])**2
            
     return result/Ns**2 - 1
     
def change_parameters(params, H, kappa, Ns) : 
      
      x = np.zeros(Ns**2, dtype = complex)
      y = np.zeros(Ns**2, dtype = complex)
      
      for i in range(Ns**2) : 
           x[i] = params[2*i] + params[2 * i + 1] * (1j)
      
      for i in range(Ns**2, 2 * Ns**2) :
           y[i-Ns**2] = params[2*i] + params[2*i + 1] * (1j)
       
      Q = params[4 * Ns**2]
      F1 = params[4 * Ns**2 + 1]
      return Energy_condensate_reduced(Q, F1, x, y, H, kappa, Ns)

def check_if_condense(H, kappa, Ns) :
    f = lambda x : change_parameters(x[0], x[1], H, kappa, Ns)
    f1 = lambda x : Energy_condensate_reduced(x[0], x[1] , np.zeros(Ns**2, dtype = complex), np.zeros(Ns**2, dtype = complex), H, kappa, Ns)    
    solucja1 = opt.minimize(f1, 5. , method = 'Nelder-Mead')    
    energy_without = solucja1.fun
    print 'energy without condensation = ', solucja1.fun, 'parameter = ', solucja1.x
#    
    
    m = 1e14
    m_params = np.zeros(4 * Ns**2 + 2)    
    for _ in range(10) : 
        params = list(np.random.rand(4*Ns**2)) + [10 * np.random.rand()] + [10 * np.random,rand()]
        solucja = opt.minimize(f, params, method = 'Nelder-Mead')
        
        print 'kappa = ', kappa, 'H = ', H, 'run = ', _ , 'energy = ', solucja.fun, 'param = ', solucja.x
        
        if solucja.fun < m :
            m = solucja.fun
            m_params = solucja.x
    
    print "final solution", m, x
    

if __name__ == "__main__":    
    
    check_if_condense(sys.argv[1], sys.argv[2], 6)