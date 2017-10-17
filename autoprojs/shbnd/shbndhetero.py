# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 15:45:14 2017

@author: leem
"""

from auto import *
import numpy as np

if __name__ == "__main__":
    my = run(e='shbnd',c='shbnd')
    

def mysave(hetero):
    print(hetero)
    timespan = hetero['PERIOD']
    A = hetero['ALPHA']
    M = hetero['M']
    N = hetero['N']
    L = hetero['LAMBDA']
    
    D = 1.0+2*A-M-N
    a0,a1 = (2+2*A-N)   / D   , (2+2*A)       / D 
    b0,b1 = (1+M)       / D   , (1+M+N)       / D 
    c0,c1 = (2+2*M)     / D   , (2+2*M+2*N)   / D 
    d0,d1 = (-2*A+2*M+N)/ D   , (-2*A+2*M+2*N)/ D 
    
    aa = a0 + a1*L
    bb = b0 + b1*L
    cc = c0 + c1*L
    dd = d0 + d1*L
    
    mystr = "ALPHA %23.15e  M %23.15e  N %23.15e LAMBDA %23.15e a %23.15e b %23.15e" %(A,M,N,L,aa,bb)
    
    np.savetxt('myshbnd.dat', np.transpose([timespan*(hetero['t']-0.5), hetero['P'],hetero['Q'],hetero['R']]))#, fmt='%23.15e',header=mystr)    
    np.savetxt('myshbnd01.dat', np.transpose([hetero['t'], hetero['P'],hetero['Q'],hetero['R']]))#, fmt='%23.15e',header=mystr) 

def changesign(data):
    for i in range(1,1001):
        data(i)['c01'] = data(i)['c01']*np.sign(data(i)['EPS0'])
        data(i)['c02'] = data(i)['c02']*np.sign(data(i)['EPS0'])