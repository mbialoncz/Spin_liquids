# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 12:47:43 2016

@author: bialonczyk
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

J=2.0

np.set_printoptions(precision=3, suppress=True)
def HamiltonianMatrixT1(n1, n2, Q, F1, F2, H, mu, kappa, Ns):
    k1 = 2*math.pi*n1/Ns
    k2 = 2*math.pi*n2/Ns
    k3=-k1-k2
    M = np.zeros((2,2),dtype=complex)
    M[0,0] = mu + H/2
    M[1,1] = mu - H/2
    M[0,1] = - J * Q[0] * (np.sin(k1)+np.sin(k2)+np.sin(k3)) * 1j
    M[1,0] =  J * Q[0] * (np.sin(k1)+np.sin(k2)+np.sin(k3)) * 1j
    
    return M

#Q[1] - X real part
#F1[1] - Z real part
def HamiltonianMatrixK01(n1, n2, Q, F1, F2, H, mu, kappa, Ns):
    k1 = math.pi*n1/Ns
    k2 = math.pi*n2/Ns
    k3 =  k1 + k2
    DH = (H/2)* np.identity(6)
    DH[3:, 3:] = -(H/2) * np.identity(3)
    DA = np.zeros((6,6), dtype = complex)
    DF = np.zeros((6,6), dtype = complex)
    DX = np.zeros((6,6), dtype = complex)
    DZ = np.zeros((6,6), dtype = complex)
    
    PA = J * Q[0]  * np.matrix([[0, - 2 * np.sin(k3) * (1j), 2 * np.sin(k1) *(1j)],
                             [-2 * np.sin(k3) * 1j, 0, 2 * np.sin(k2) * 1j],
                              [2 * np.sin(k1) * 1j, 2 * np.sin(k2) * 1j, 0]])
    
    PF1 = J * F1[0] * np.matrix([[0, 2* np.cos(k3), 2 * np.cos(k1)],
                               [2 * np.cos(k3), 0, 2 * np.cos(k2)],
                                [2 * np.cos(k1), 2 * np.cos(k2), 0]])
                                
    PF2 = 1j * J * F2 * np.matrix([[0, 2 * np.cos(k3), -2 * np.cos(k1)], 
                                   [-2 * np.cos(k3), 0, 2 * np.cos(k2)], 
                                     [2 * np.cos(k1), -2 * np.cos(k2), 0]])
                            
    PX = J * Q[1] * np.matrix([[0, - 2 * np.sin(k3) * (1j),  -2 * np.sin(k1) *(1j)],
                             [2 * np.sin(k3) * 1j, 0, 2 * np.sin(k2) * 1j],
                              [2 * np.sin(k1) * 1j, - 2 * np.sin(k2) * 1j, 0]])
                              
                             
    PZ1 = J * F1[1] * np.matrix([[0, 2 * np.cos(k3), -2 * np.cos(k1)], 
                                   [2 * np.cos(k3), 0, 2 * np.cos(k2)], 
                                     [-2 * np.cos(k1), 2 * np.cos(k2), 0]])
    
    DA[:3, 3:] = PA
    DA[3:, :3] = PA.H
    
    DF[:3, :3] = PF1 + PF2
    DF[3:, 3:] = PF1 + PF2
    
    DX[:3, 3:] = PX
    DX[3:, :3] = PX.H
    
    DZ[:3, 3:] = PZ1
    DZ[3:, :3] = -PZ1
    
    
    Dmu = mu * np.identity(6)
    
    return - 1/2 *  DA - DH + Dmu  
    
    

def HamiltonianMatrixK02(n1, n2, Q, F1, F2, H, mu, kappa, Ns):
    k1 = math.pi*n1/Ns
    k2 = math.pi*n2/Ns
    k3 = k1 + k2 
    DH = (H/2) * np.identity(6)
    DH[3:, 3:] = - (H/2) * np.identity(3)
    DA = np.zeros((6,6), dtype = complex)
    DF = np.zeros((6,6), dtype = complex)    
    
    
    PA = J *  Q[0] * np.matrix([[0, 2 * np.cos(k3), 2 * np.cos(k1)], 
                [-2 * np.cos(k3), 0, 2 * np.cos(k2)], 
                    [-2 * np.cos(k1), -2 * np.cos(k2), 0]]) 

    PF1 = J * F1[0] * np.matrix([[0, 2 * np.cos(k3), -2 * 1j * np.cos(k2)], 
                  [2 * np.cos(k3), 0, 2 * np.cos(k2)], 
                    [2 * 1j * np.cos(k2), 2 * np.cos(k2), 0]]) ;

    PF2 = J * 1j * F2 * np.matrix([[0, 2 * np.cos(k3), 2 * 1j * np.cos(k2)], 
                                    [-2 * np.cos(k3), 0, 2 * np.cos(k2)], 
                                       [2 * 1j * np.cos(k2), -2 * np.cos(k2), 0]]) 
                                       
    DA[:3, 3:] = PA
    DA[3:, :3] = PA.H
    
    DF[:3, :3] = PF1 + PF2
    DF[3:, 3:] = PF1 + PF2
    
    Dmu = mu * np.identity(6)
    
    return 1/4 * (DF - DA) - DH + Dmu
    
def HamiltonianMatrixSa1(n1, n2, Q, F1, F2, H, mu, kappa, Ns) :
    k1 =  math.pi*n1/Ns
    k2 =  math.pi*n2/Ns
    k3 = - k1 - k2 
    DH = H * np.identity(6)
    DH[3:, 3:] = - H * np.identity(3)
    DA = np.zeros((6,6), dtype = complex)
    
#    PA1 = J * Q * np.matrix([[0, np.exp(-k1* 1j)+ np.exp(k1 * 1j),-np.exp(-k3* 1j)-np.exp(k3 * 1j)],
#                              [-np.exp(-k1* 1j) - np.exp(k1 * 1j), 0, np.exp(-k2* 1j)+ np.exp(k2 * 1j) ],
#                                      [np.exp(-k3 * 1j)+ np.exp(k3 * 1j), -np.exp(-k2 * 1j)- np.exp(k2 * 1j), 0]])
#    
#    PA = J * Q * np.matrix([[0, np.cos(k1),-np.cos(k3)],
#                                [-np.cos(k1), 0, np.cos(k2) ],
#                                     [np.cos(k3), -np.cos(k2), 0]])
                                     
    PA2 = J * np.matrix([[0, Q[0] *np.exp(-k1* 1j)+ Q[1] * np.exp(k1 * 1j),- Q[1] * np.exp(-k3* 1j) - Q[0] * np.exp(k3 * 1j)],
                              [-Q[1] * np.exp(-k1* 1j) - Q[0] * np.exp(k1 * 1j), 0, Q[0] * np.exp(-k2* 1j)+ Q[1] * np.exp(k2 * 1j)],
                                      [Q[0] * np.exp(-k3 * 1j)+ Q[1] * np.exp(k3 * 1j), - Q[1] * np.exp(-k2 * 1j)- Q[0] * np.exp(k2 * 1j), 0]])
    
    DA[:3, 3:] = PA2
    DA[3:, :3] = PA2.H
    
    Dmu = mu * np.identity(6)
#    print PA
#    print PA1
#    print 1/2 * DA + Dmu - DH
    return 1/2 * DA + Dmu - DH

def HamiltonianMatrixSa2(n1, n2, Q, F1, F2, H, mu, kappa, Ns) :
    k1 = 2 * math.pi*n1/Ns
    k2 = 2 * math.pi*n2/Ns
    k3 = - k1- k2 
    
    DH = H * np.identity(6)
    DH[3:, 3:] = - H * np.identity(3)    
    
    
    DA = np.zeros((6,6), dtype = complex)
    PA = J * Q * np.matrix([[0, np.exp(-k1* 1j) - np.exp(k1 * 1j),-np.exp(-k3* 1j) + np.exp(k3 * 1j)],
                              [np.exp(-k1* 1j) - np.exp(k1 * 1j), 0, np.exp(-k2* 1j) - np.exp(k2 * 1j) ],
                                      [np.exp(-k3* 1j)- np.exp(k3 * 1j), np.exp(-k2 * 1j)- np.exp(k2 * 1j), 0]])    
    
    
    DA[:3, 3:] = PA
    DA[3:, :3] = PA.H
    
    Dmu = mu * np.identity(6)
    
    return 1/2 * DA + Dmu - DH
    
def HamiltonianMatrixKpi1(n1, n2, Q, F1, F2, H, mu, kappa, Ns):
    k1 = 2 * math.pi*n1/Ns
    k2 = 2 * math.pi*n2/Ns
    k3 = k1 + k2 
    DA = np.zeros((12,12), dtype = complex)
    
    PA = J * Q * np.matrix([[0, -np.exp(k3 * 1j), 2*np.sin(k1) * 1j, 0, np.exp(-k3 * 1j),0],
                             [np.exp(-k3 * 1j), 0, -np.exp(-k2 * 1j), -np.exp(k3 * 1j),0, np.exp(k2 * 1j)],
                              [2 * np.sin(k1) * 1j, np.exp(k2 * 1j), 0, 0, np.exp(k2 * 1j), 0],
                               [0, np.exp(-k3 * 1j), 0, 0, np.exp(k3 * 1j), -2*np.cos(k1)],
                                [-np.exp(k3 * 1j), 0 , np.exp(k2 * 1j), - np.exp(-k3 * 1j), 0, -np.exp(-k2 * 1j)],
                                 [0, np.exp(-k2 * 1j), 0, 2 * np.cos(k1), -np.exp(k2 * 1j), 0]])
    DA[:6, 6:] = PA
    DA[6:, :6] = PA.H    
    
    Dmu = mu * np.identity(12)
    DH = (H/2) * np.identity(12)
    DH[6:,6:] = -H * np.identity(6)
    
    return -1/2 * DA - DH + Dmu



def HamiltonianMatrixKpi2(n1, n2, Q, F1, F2, H, mu, kappa, Ns):
    
    k1 = 2 * math.pi*n1/Ns
    k2 = 2 * math.pi*n2/Ns
    k3 = k1 + k2 
    DA = np.zeros((12,12), dtype = complex)
    
    PA = J * Q * np.matrix([[0, np.exp(k3 * 1j), -2*np.cos(k1) * 1j, 0, np.exp(-k3 * 1j), 0],
                             [-np.exp(-k3 * 1j), 0, np.exp(-k2 * 1j), -np.exp(k3 * 1j),0, np.exp(k2 * 1j)],
                              [2 * np.cos(k1) * 1j, -np.exp(k2 * 1j), 0, 0, -np.exp(k2 * 1j), 0],
                               [0, np.exp(-k3 * 1j), 0, 0, -np.exp(k3 * 1j), -2*np.sin(k1) * 1j],
                                [-np.exp(k3 * 1j), 0 , np.exp(k2 * 1j), - np.exp(-k3 * 1j), 0, np.exp(-k2 * 1j)],
                                 [0, -np.exp(-k2 * 1j), 0, -2 * np.sin(k1)*1j, -np.exp(k2 * 1j), 0]])
    DA[:6, 6:] = PA
    DA[6:, :6] = PA.H    
    
    Dmu = mu * np.identity(12)
    DH = (H/2) * np.identity(12)
    DH[6:,6:] = -H * np.identity(6)
    
    return -1/2 * DA - DH + Dmu

    
def HamiltonianMatrix(n1, n2, Q, F1, F2, H, mu, kappa, Ns, ansatz) :
    
    if ansatz == 'T1':
        return [HamiltonianMatrixT1(n1, n2, Q, F1, F2, H, mu, kappa, Ns),2]
    elif ansatz == 'K01' :
        return [HamiltonianMatrixK01(n1, n2, Q, F1, F2, H, mu, kappa, Ns),6]
        
    elif ansatz == 'K02' :
        return [HamiltonianMatrixK02(n1, n2, Q, F1, F2, H, mu, kappa, Ns),6]
    
    elif ansatz == 'Sa1' :
        return [HamiltonianMatrixSa1(n1, n2, Q, F1, F2, H, mu, kappa, Ns),6]
    elif ansatz == 'Sa2' :
        return [HamiltonianMatrixSa1(n1, n2, Q, F1, F2, H, mu, kappa, Ns),6]
    
    
    elif ansatz == 'Kpi1' :
        return [HamiltonianMatrixKpi1(n1, n2, Q, F1, F2, H, mu, kappa, Ns),12]
    
    elif ansatz == 'Kpi2' :
        return [HamiltonianMatrixKpi2(n1, n2, Q, F1, F2, H, mu, kappa, Ns),12]
        

                  
def Energy(Q, F1, F2, H, mu, kappa, Ns, ansatz) :
    
    #print Q
    result = 0

    for n1 in range(Ns) :
        for n2 in range(Ns) :

            M, dim = HamiltonianMatrix(n1, n2, Q, F1, F2, H, mu, kappa, Ns, ansatz) 
        
            B = np.identity(dim)
            B[dim/2:dim, dim/2:dim] = -np.identity(dim/2)
        
            eig = np.absolute(np.real(lin.eigvals(np.dot(B,M))))
        
            result += sum(eig)/2
    
    if ansatz == 'T1' :
        return result/(Ns**2) +J/2 * 3 *  Q[0]**2  - mu*(1. + kappa) - H/2
    elif ansatz == 'Sa1' or ansatz == 'Sa2' :
        return result/(3 * Ns**2) + J * Q**2 - mu*(1.+kappa) - H/2
    elif ansatz == 'K01' or ansatz == 'K02' : 
        return result/(3 * Ns**2) + J * Q[0]**2 - mu * (1.+kappa) - H/2
        

def dispersion(k1, k2, Q, F1, F2, H, mu, kappa, Ns, ansatz) :
    
    n1 = Ns * k1/ (2*math.pi)
    n2 = Ns * k2/ (2*math.pi)
    M, dim = HamiltonianMatrix(n1, n2, Q, F1, F2, H, mu, kappa, Ns, ansatz) 

    B = np.identity(dim)
    B[dim/2:dim, dim/2:dim] = -np.identity(dim/2)

    eig = np.absolute(np.real(lin.eigvals(np.dot(B,M))))

    return min(eig)
            
    
    

def BogolubovTransformation(n1, n2, Q, F1, F2, H, mu, kappa, Ns, ansatz) :
    
     M, dim = HamiltonianMatrix(n1, n2, Q, F1, F2, H, mu, kappa, Ns, ansatz)
    
     B = np.matrix(np.identity(dim))
     B[dim/2:dim, dim/2:dim] = -np.identity(dim/2)    
     
     try :
        u,w = lin.eig(np.dot(B,M))
     except np.linalg.LinAlgError as err :
         print M
         print n1, n2, Q, F1, F2, H, mu
         raise lin.LinAlgError
     w=np.matrix(w)     
     
     temp = np.real(np.dot(np.dot(w.H,B),w))
    # print np.diagonal(temp)
     order = np.argsort(np.diagonal(temp))[::-1]     
     
     w = (w.transpose()[order]).transpose()

     k1 = [(np.square(np.absolute(w[:dim//2,i]))).sum()-(np.square(np.absolute(w[dim//2: ,i]))).sum() for i in range(dim//2)]
     k2 = [-(np.square(np.absolute(w[:dim//2,i]))).sum()+(np.square(np.absolute(w[dim//2 : ,i]))).sum() for i in range(dim//2,dim)]
    
     k = k1 + k2 
     

     U = np.matrix(np.zeros((dim, dim),dtype=complex))  
    
     
     for i in range(dim//2) :
         U[:,i] = w[:,i]/(np.sqrt(k[i]))
         
     for i in range(dim//2,dim):
         U[:,i] = w[:,i]/np.sqrt(k[i])
         
     
#     print np.absolute(U[:dim/2,dim/2:]-U[dim/2:,:dim/2]
#     print np.diagonal(np.dot(np.dot(U.H,M),U))
#     print np.diagonal(np.dot(np.dot(lin.inv(U),np.dot(B,M)),U))
     l1 = [(np.square(np.absolute(U[:dim//2,i]))).sum()-(np.square(np.absolute(U[dim//2: ,i]))).sum() for i in range(dim//2)]
     l2 = [-(np.square(np.absolute(U[:dim//2,i]))).sum()+(np.square(np.absolute(U[dim//2 : ,i]))).sum() for i in range(dim//2,dim)]
#     print l1+l2
     return U
     
def ParticleNumber(Q, F1, F2, H, mu, kappa, Ns, ansatz) :
    
    if ansatz == 'T1' :
        dim, factor = 2, 1 
    elif ansatz == 'K01' or ansatz == 'K02' or ansatz == 'Sa1' or ansatz == 'Sa2' :
        dim, factor = 6, 3
    else :
        dim, factor = 12, 3
    
    result = 0
    for n1 in range(Ns) :
        for n2 in range(Ns) :
            
            U = BogolubovTransformation(n1, n2, Q, F1, F2, H, mu, kappa, Ns, ansatz)
            
            U12 = np.absolute(U[dim/2:,:dim/2])**2
            U21 = np.absolute(U[:dim/2,dim/2:])**2

            result += U12.sum()+U21.sum()
     
    return np.real(result/(factor * Ns**2))
    
def AverageMagnetization(Q, F1, F2, H, mu, kappa, Ns, ansatz) :
    
    if ansatz == 'T1' :
        dim, factor = 2,1
    elif ansatz in set(['K01','K02','Sa1','Sa2']) :
        dim, factor = 6, 3
    else :
        dim = 12
    
    result = 0
    for n1 in range(Ns) :
        for n2 in range(Ns) :
            
            U = BogolubovTransformation(n1, n2, Q, F1, F2, H, mu, kappa, Ns, ansatz)
            
            U12 = np.absolute(U[dim/2:,:dim/2])**2
            U21 = np.absolute(U[:dim/2,dim/2:])**2
            
            result += U12.sum()-U21.sum()
            
    return result/(2*factor * Ns**2)  + H/2   
    

def Bis(a, b, Q, F1, F2, H, kappa, Ns, ansatz) :
    
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

    while (ParticleNumber(Q, F1, F2, H, b, kappa,Ns,ansatz) == np.inf or
            np.isnan(ParticleNumber(Q, F1, F2, H, b, kappa,Ns,ansatz)) or
              ParticleNumber(Q, F1, F2, H, b, kappa,Ns,ansatz) >= kappa):
                  b = 2 * b
    
    while (ParticleNumber(Q, F1, F2, H, a, kappa,Ns,ansatz) == np.inf or 
            np.isnan(ParticleNumber(Q, F1, F2, H, a, kappa,Ns,ansatz))  or 
                np.absolute(ParticleNumber(Q, F1, F2, H, a, kappa,Ns,ansatz) -  
                      ParticleNumber(Q,F1,F2,H,b, kappa,Ns,ansatz)) > 0.00001) :
        
        c = (a+b)/2
        it += 1
        
        if it > 30:
            return c
        
        #print Q, c
        
        if (np.sign(ParticleNumber(Q, F1, F2, H, c, kappa, Ns, ansatz) - kappa) == 
                np.sign(ParticleNumber(Q, F1, F2, H, b, kappa, Ns, ansatz)-kappa)) :
            b = c
        else : 
            a = c
            
    return c
    
def Metropolis() :
    pass
    
    
    
def Energia(a, b, Q, F1, F2, H,  kappa, Ns, ansatz) :
    try :
        
         mu = Bis(a, b, Q, F1, F2, H, kappa, Ns, ansatz)
#         print Q, mu
    except np.linalg.LinAlgError :
        return 1e14
    
    
    return Energy(Q, F1, F2, H, mu, kappa, Ns, ansatz)

def spectral_gap(Q, F1, F2, H,  kappa, mu,  Ns, ansatz) :
    
    energies_list = []
        
    for n1 in range(Ns) :
        for n2 in range(Ns) :
            
            M, dim = HamiltonianMatrix(n1, n2, Q, F1, F2, H, mu, kappa, Ns, ansatz)
            B = np.identity(dim)
            B[dim/2:dim, dim/2:dim] = -np.identity(dim/2)
        
            eig = np.absolute(np.real(lin.eigvals(np.dot(B,M))))
            energies_list = np.concatenate((eig, energies_list), 0)
            
            
    return min(energies_list)
    
def lower_branch(k1, k2, Q, F1, F2, H, mu, kappa, Ns, ansatz):
    M, dim = HamiltonianMatrix(n1, n2, Q, F1, F2, H, mu, kappa, Ns, ansatz)
    B = np.identity(dim)
    B[dim/2:dim, dim/2:dim] = -np.identity(dim/2)
    
    eig = np.absolute(np.real(lin.eigvals(np.dot(B,M))))
    return min(eig)
    
    
    
    
    
    
    
    
    
    
    
    
    
    

