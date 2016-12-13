# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 10:53:39 2016

@author: mbialoncz
"""

Ns = 12

with open('result_K02_12', 'r') as res :
    results = []
    for line in res : 
        print line
        results.append([float(x) for x in line.split()])
        
kappa = results[1][0]
Q_minimal = results[1][1]
mu_minimal = results[1][2]

print modes_of_dispersion([Q_minimal,0], [0,0], 0, 0, mu_minimal, kappa, Ns, 'K02')