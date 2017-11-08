# -*- coding: utf-8 -*-
"""
Created on Sun Oct  8 12:05:17 2017

@author: leem
"""

import matplotlib.pyplot as myplot
import scipy as sp
import numpy as np
#from scipy.integrate import odeint

# this data was saved from auto. It was input in the Fortran code, 
# and safer not computing in python but importing pre-existing values.
mypar = np.load('mypar.npy')

# 4 eigenvectors at M0
X01 = np.zeros(4)
X02 = np.zeros(4)
X03 = np.zeros(4)
X04 = np.zeros(4)

# 4 eigenvectors at M1
X11 = np.zeros(4)
X12 = np.zeros(4)
X14 = np.zeros(4)
X13 = np.zeros(4)

# those in W^u(M0), W^s(M1)
X01 = mypar[21:25]
X02 = mypar[25:29]
X03 = mypar[29:33]
X11 = mypar[33:37]
X12 = mypar[37:41]
X14 = mypar[41:45]

# those in W^s(M0), W^u(M1). not really necessary.
X13 = mypar[45:49]
X04 = mypar[49:53]  

TMAX = 10.0
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
#E0= C1
#F0= 1.0
E0 =  C1/(C1+np.exp(TMAX))
F0 = 1.0/(C1+np.exp(TMAX))



def myp(t):
    return 0.0
    
def myq(t):
    C1=1.0
    return 1.0/(1.0+C1*np.exp(-t)) 

def myr(t):
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
# for given p(t),q(t),r(t), compute the r-h-s f(t,x) = f(t,p(t),q(t),r(t),s)    
    P = myp(t)
    Q = myq(t)
    R = myr(t)
    return ( S * ( (R-aa)*(A-M-N)/(L*(1.0+A)) + L*P*R + Q -   ( R*( S-(1.0+M+N)/(1.0+A) ) + N/(1.0+A) )  /L  ) )
    

# this module integrates the s data in [-TMAX,TMAX]  
 
# compute initial point

MAT = np.zeros((3,3))
MAT = np.transpose([X01[0:3],X02[0:3],X03[0:3]])
myb = np.transpose([myp(-TMAX)-0.0, myq(-TMAX)-0.0,myr(-TMAX)-R0])
myc0 = np.linalg.solve(MAT,myb)
eps0 = np.linalg.norm(myc0)
myc0 = myc0/eps0
u0 = [0.0,0.0,R0,S0] + eps0*( myc0[0]*X01 +  myc0[1]*X02 +myc0[2]*X03  )    

# integrate s   
myint = sp.integrate.ode(myf).set_integrator('dopri5')    
myint.set_initial_value(u0[3], -TMAX)

t = np.linspace(-TMAX,TMAX,200)  
P=np.zeros(t.size)  
Q=np.zeros(t.size)
R=np.zeros(t.size)
S=np.zeros(t.size)
k=t[1]-t[0]

P[0],Q[0],R[0],S[0]=u0
for i in range(1,200):
    P[i],Q[i],R[i] = myp(t[i]), myq(t[i]), myr(t[i])
    S[i]=myint.integrate(t[i])    
    if not myint.successful(): break

# compute final point and 

MAT = np.zeros((4,4))
MAT = np.transpose([X11,X12,X13,X14])
myb = np.transpose([P[-1]-0.0, Q[-1]-1.0,R[-1]-R1, S[-1]-S1])
myc1 = np.linalg.solve(MAT,myb)
eps1 = np.linalg.norm(myc1)
myc1 = myc1/eps1
print(myc1)
print(eps1)

np.savetxt('starting_solution_scaled.dat', np.transpose([t/20.0+0.5, P, Q , R, S]), fmt='%23.15E')       
    
if __name__=="__main__":
#    t = np.linspace(-10,10,100)
#    R = np.zeros(t.size)
#    for i in range(100):
#        R[i]=myr(t[i])
    myplot.plot(t,R)
    S_fig=myplot.figure()
    S_figp = S_fig.add_subplot(111)
    S_figp.plot(t,S)
        
        
def testing_formulas():        
    QA = 1.0
    QB = -( (1-S0)/L - N/(L*R0) )*R0/N + S0*R0/L
    QC = -( S0*R0*R0/(L*N) )*( (1-S0)/L - N/(L*R0) ) - (S0*R0/N)*((1-S0)/L)*(A*R0/L)
    
    mu01 = 2.0
    mu02 = 1.0
    mu03 = (  -QB + np.sqrt(QB**2-4*QA*QC)  )/(2*QA)
    mu04 = (  -QB - np.sqrt(QB**2-4*QA*QC)  )/(2*QA)
    print(mu03+(M+N)*aa/(N*L))
    
    
    QA = 1.0
    QB = -( (1.0-S1)/L - N/(L*R1) )*R1/N + S1*R1/L
    QC = -( S1*R1*R1/(L*N) )*( (1-S1)/L - N/(L*R1) ) - (S1*R1/N)*((1-S1)/L)*(A*R1/L)
    mu11 = -(1.0+M+N)/(A-M-N)
    mu12 = -1.0
    mu13 = (  -QB + np.sqrt(QB**2-4*QA*QC)  )/(2*QA)
    mu14 = (  -QB - np.sqrt(QB**2-4*QA*QC)  )/(2*QA)
    print(mu13+(M+N)*R1/(N*L))
    
    MAT0 = np.zeros((4,4))
    MAT1 = np.zeros((4,4))
    
    MAT0 = [ [2, 0, 0, 0], [bb*R0, 1, 0, 0], [(R0/N)*L*R0, R0/N, (R0/N)*( (A-M-N)/(L*(1.0+A)) - N*A/( L*(1.0+A)*R0 )  ), (R0/N)*A*R0/L], [S0*L*R0, S0, S0*( (A-M-N)/(L*(1.0+A)) + N/(L*(1.0+A)*R0)   ), S0*(-R0)/L]]
    MAT1 = [[-(1+M+N)/(A-M-N), 0, 0, 0], [(bb-L)*R1, -1, 0, 0], [(R1/N)*L*R1, R1/N, (R1/N)*( (A-M-N)/(L*(1.0+A)) - N*A/(L*(1.0+A)*R1)   ), (R1/N)*A*R1/L], [S1*L*R1, S1, S1*( (A-M-N)/(L*(1.0+A)) + N/(L*(1.0+A)*R1)   ), S1*(-R1)/L]]       

