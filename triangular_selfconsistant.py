
from __future__ import division
import numpy as np
import scipy.linalg as lin
import scipy.optimize as opt
from scipy.integrate import quad, dblquad
import math
import random
import matplotlib.pyplot as plt
import sys
import random as rand

from mpl_toolkits.mplot3d import Axes3D
#from mean_field import *
from pylab import *
from mean_field import *

J=1.
B=0.05
r = 3. * sqrt(3) / 2. 

N=20



def Ek(k, mu, delta) :
#    print (mu ** 2 -  delta**2 * J**2 * (np.sin(k[0])+np.sin(k[1]))**2)**0.5
    k1 = k[0]
    k2 = k[1]
    k3 = k1 + k2
    return (mu ** 2 -  4 * delta**2 * J**2 * (np.sin(k[0])+np.sin(k[1]) - np.sin(k3))**2)**0.5
    
k_val = np.linspace(-np.pi,np.pi,N)
k_vals = []
for k1 in k_val :
    for k2 in k_val :
        k_vals.append([k1,k2])
        

"""
x = [lambda, delta, m]
"""

def self1(x) :
    R = 0
    for k in k_vals :
        R += x[0]/Ek(k,x[0],x[1])
    
    return R/N**2 + 2 * x[0] * x[2]/(1/2 * B - 2 * x[2] * J)

def self2(x) : 
    R = 0
    for k in k_vals : 
        R +=  (np.sin(k[0]) + np.sin(k[1])-sin(k[0]+k[1]))**2/Ek(k,x[0],x[1])
    
    return R* J /(3 * N**2) + 2*r**2/3. * J * x[2]/(1/2 * B - 2 * x[2] * J)
    
def self3(x) : 
    return Ek([2*np.pi/3, 2 * np.pi/3],x[0], x[1]) - (1/2 * B - 2 * x[2] * J)
    

f = lambda x : [self1(x)-2, self2(x)-1, self3(x)] 

B= 7.5
x0 = [1/2 + r**2/3 ,0.0, 0.5]
print self1(x0), self2(x0), self3(x0)

B_values = []
m_values = []
delta_values = []

x0 = [2.48, 0.05, 0.51]
for B in np.arange(7.5,0,-0.1) :
    sol = self_consistant(f, x0)
    B_values.append(B)
    m_values.append(sol[2])
    delta_values.append(sol[1])
    print B, sol[2] 
    x0 = sol

plt.plot(B_values, m_values)
plt.show()
    
    
    
    
    
    
    

        
