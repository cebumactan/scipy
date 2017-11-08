# -*- coding: utf-8 -*-
"""
Created on Sun Oct  8 12:05:17 2017

@author: leem
"""

#from numpy import *
#from scipy import *
import matplotlib.pyplot as myplot
import scipy as sp
import numpy as np
#import scipy.io
#from scipy.integrate import odeint

A = 0
M = -0.5
K = 10
N = 1.0/K
LMAX = 2*(A-M-N)*(1+M)/((1+M+N)**2)
L = 0.5*LMAX

D = 1.0+2*A-M-N
a0,a1 = (2+2*A-N)   / D   , (2+2*A)       / D 
b0,b1 = (1+M)       / D   , (1+M+N)       / D 
c0,c1 = (2+2*M)     / D   , (2+2*M+2*N)   / D 
d0,d1 = (-2*A+2*M+N)/ D   , (-2*A+2*M+2*N)/ D 

aa = a0 + a1*L
bb = b0 + b1*L
cc = c0 + c1*L
dd = d0 + d1*L

R0 = aa
R1 = R0 - (1.0+A)*L/(A-M-N)
S0 = (1.0+M+N)/(1.0+A) - N/(R0*(1.0+A))
S1 = (1.0+M+N)/(1.0+A) - N/(R1*(1.0+A))

C1 = 1.0
#    E0 =  C1/(C1+np.exp(TMAX))
#    F0 = 1.0/(C1+np.exp(TMAX))
E0= C1
F0= 1.0


def myp(t):
    return 0.0
    
def myq(t):
    C1=1.0
    return 1.0/(1.0+C1*np.exp(-t)) 

def myr(t):
#    RM = R1 - 1e-4
    W = -(M+N)*aa/L
    NUMERATOR = aa*(E0+F0*np.exp(t))**K

    DENOM1 = 0
    for j in range(0,K+1):
        DENOM1 = DENOM1 + sp.misc.comb(K,j)*(E0**(K-j))*(F0**j)*(np.exp(j*t)) / (-K*W+j)


    DENOM1 = (-K*W) * DENOM1

#    for j in range(0,K+1):
#        DENOM2 = DENOM2 + comb(K,j)*(E0**(K-j))*(F0**j)*(np.exp(j*TMAX)) / (-K*W+j)
#
#    DENOM2 = ( 1 + (K*W*RM/aa)*DENOM2 ) * np.exp(-K*W*(t-M))
#
#    R = NUMERATOR / (DENOM1 + 0*DENOM2)    
    
    R = NUMERATOR / DENOM1
    return R
    
def myf(t,S):
#    A = 0
#    M = -0.5
#    K = 10
#    N = 1.0/K
#    LMAX = 2*(A-M-N)*(1+M)/((1+M+N)**2)
#    L = 0.5*LMAX
#    
#    D = 1.0+2*A-M-N
#    a0,a1 = (2+2*A-N)   / D   , (2+2*A)       / D 
#    b0,b1 = (1+M)       / D   , (1+M+N)       / D 
#    c0,c1 = (2+2*M)     / D   , (2+2*M+2*N)   / D 
#    d0,d1 = (-2*A+2*M+N)/ D   , (-2*A+2*M+2*N)/ D 
#    
#    aa = a0 + a1*L
#    bb = b0 + b1*L
#    cc = c0 + c1*L
#    dd = d0 + d1*L
#    
#    R0 = aa
#    R1 = R0 - (1.0+A)*L/(A-M-N)
#    S0 = (1.0+M+N)/(1.0+A) - N/(R0*(1.0+A))
#    S1 = (1.0+M+N)/(1.0+A) - N/(R1*(1.0+A))
    
    P = myp(t)
    Q = myq(t)
    R = myr(t)
    return ( S * ( (R-aa)*(A-M-N)/(L*(1.0+A)) + L*P*R + Q -   ( R*( S-(1.0+M+N)/(1.0+A) ) + N/(1.0+A) )  /L  ) )
    
if __name__=="__main__":
    t = np.linspace(-10,10,100)
    R = np.zeros(t.size)
    for i in range(100):
        R[i]=myr(t[i])
    myplot.plot(t,R)
        
    

