# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 10:27:34 2017

@author: mbialoncz
"""

# -*- coding: utf-8 -*-


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
import sys

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
    f = lambda x : change_parameters(x, H, kappa, Ns)
#    f1 = lambda x : Energy_condensate_reduced(x[0], 0 , np.zeros(Ns**2, dtype = complex), np.zeros(Ns**2, dtype = complex), H, kappa, Ns)    
#    solucja1 = opt.minimize(f1, 5. , method = 'Nelder-Mead')    
#    energy_without = solucja1.fun
#    print 'energy without condensation = ', solucja1.fun, 'parameter = ', solucja1.x
#    
    
    m = 1e14
    m_params = np.zeros(4 * Ns**2 + 2)    
    for _ in range(10) : 
        params = list(np.random.rand(4*Ns**2)) + [10 * np.random.rand()] + [0.]
        solucja = opt.minimize(f, params, method = 'Nelder-Mead')
        
        print 'kappa = ', kappa, 'H = ', H, 'run = ', _ , 'energy = ', solucja.fun, 'param = ', solucja.x
        
        if solucja.fun < m :
            m = solucja.fun
            m_params = solucja.x
    
    print "final solution", m, x
    
if __name__ == "__main__":    
    
    check_if_condense(sys.argv[1], sys.argv[2], 6)
        
        
    
    
#x0 = np.zeros(N**2, dtype = complex)
#y0 = np.zeros(N**2, dtype = complex)
#x1 = np.zeros(N**2, dtype = complex)
#y1 = np.zeros(N**2, dtype = complex)

#for n1 in range(N) :
#    for n2 in range(N) :
#        x0[N * n1 + n2] = 0
#        y0[N * n1 + n2] = 0

#x0[N**2/3 + N/3] = np.random.rand() + np.random.rand()*1j
#y0[N**2/3 + N/3] = np.random.rand() + np.random.rand()*1j
#Q0 = 0.2
#F10 = 0.
#mu0 = 2.5
#x0[N**2/3 + N/3] = 1. 
#y0[N**2/3 + N/3] = 1j
#x0[2/3 * N**2 + 2/3 * N] = 1
#y0[2/3 * N**2 + 2/3 * N] = -1j 

#print find_minimum(Q0, 0.49, mu0 ,kappa, N)
#print Average_number(Q0, F10, x1, y1, 0., mu0, kappa, N)
#print Average_number1(Q0, F10, x0, y0, 0., mu0, )
#f1 = lambda params1 : Energy_condensate([params1[0], 0, params1[1],params1[2],params1[3],params1[4],params1[5],params1[6],params1[7],params1[8]],0,kappa,N)
#params0 = [0.2, 1., 0, 0, 0, 0, 0, 0, 1, 0]
#
#for kappa in np.arange(0.2, 1.4, 0.2) :
#    m=1000000000
#    po = np.zeros(9)
#    for _ in range(20) :
#    
#           params01 = [5.] + list(np.random.rand(8))
#           print params01
#           solucja = opt.minimize(f1, params01, method = 'Nelder-Mead')
#           print 'sol' , _, solucja.x, solucja.fun
#           if solucja.fun <= m :
#               po = solucja.x
#               m = solucja.fun
#    print kappa, m, po
#
#print 'result', m
#
#f = lambda params : Energy_condensate(params, 0., kappa, N)

#
#solucja = opt.minimize(f1, params01, method = 'Nelder-Mead')
#print 'sol' , solucja.x
#print Average_number(Q0, F10, x0, y0, 0., a, kappa, N) 
#
#mu_values = np.arange(1.2,1.5,0.05)
#print Average_number(Q0, F10,x0,y0,0,1.0311, kappa, N)
#print Bis_condensate(0, 5, Q0, 0, x0, y0, 0, kappa, N)
#
#Q_values = np.arange(-6,6, 0.1)
#F_values = np.arange(-0.5,0.5, 0.1)
#particle_values =  [Average_number(Q0, F10, x0, y0, 0., mu, kappa, N) for mu in mu_values]
#particle_values1 = [Average_number1(Q0, F10, x0, y0, 0., mu, kappa, N) for mu in mu_values]
#Energy_values = [np.abs(Energy_condensate_reduced(Q, 0, x0, y0, 0.0, kappa, N)) for Q in Q_values]
##print particle_values
#print Energy_values
#mu_v = [Bis_condensate(0,5,Q,0,x0,y0,0,kappa, N) for Q in Q_values]
#print Bis_condensate(0,5,Q,0,x0,y0,0,kappa, N)
#print 'av', Average_number1(Q, 0.0,x0,y0, 0, Bis_condensate(0,5,Q,0,x0,y0,0,kappa, N), kappa,N)
#Energy_values = [np.abs(f1([Q])) for Q in Q_values]
#plt.plot(Q_values, Energy_values)
#plt.savefig('dupa1.jpg') 
     

     
     
           
    
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
 
            
    