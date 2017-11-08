#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 15:38:02 2017, KAUST

This simple code does the following : 
    
given the values of (p,q,r) at discrete "times" t 

it integrates numerically the 4th equation is the (p,q,r,s) developed by Min-Gi

 
@author: Thodoros, Min-Gi
"""


from math import *
import numpy as np

#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

#Machine epsilon
eps = 1.0;                         
while 1 + eps != 1: eps = eps/2.0

atol = 1.e-9
rtol = 1.e-3

## Constants for the Runge-Kutta 4-5 method
a21=0.25;                                                                                    c2=0.25;
a31=3.0/32.0;      a32=9.0/32.0;                                                             c3=3.0/8.0;
a41=1932.0/2197.0; a42=-7200.0/2197.0; a43=7296.0/2197.0;                                    c4=12.0/13.0;
a51=439.0/216.0;   a52=-8.0;           a53=3680.0/513.0;    a54=-845.0/4104.0;               c5=1.0; 
a61=-8.0/27.0;     a62=2.0;            a63=-3544.0/2565.0;  a64=1859.0/4104.0; a65=-11.0/40; c6=0.5;
#--------------------------------------------------------------------------
b1=25.0/216.0;     b3=1408.0/2565.0;   b4=2197.0/4101.0;    b5=-0.2;  # low order method
bb1=16.0/135.0;    bb3=6656.0/12825.0; bb4=28561.0/56430.0; bb5=-9.0/50.0; bb6=2.0/55.0;   #high order method
#--------------------------------------------------------------------------

#now define the r.h.s for the s-equation 
def F(s,params):
    p=params[0];   q=params[1]; r=params[2]
    alf=params[3]; m=params[4]; n=params[5]
    lam=params[6]; a=params[7];
    
    t1=(alf-m-n)*(r-a)/(lam*(1.0+alf))
    t2=lam * p * r 
    t3=r*(s - (1.0+m+n)/(1.0+alf))/lam
    t4=n/(lam*(1.0+alf))
    
    return s*(t1+t2+q-t3-t4)*20.0;

# Copy sign of b to a (assuming a > 0).
def sign(a,b):                  
  if b < 0: return -a
  else: return a  


  
  
def rkf45(f,y,a,h,da, params):
# f:      subroutine to compute y_prime: f(x,y,yp) where x = independent
#         variable, y = y[...] is the vector of dependent variables, and
#         yp = yp[...] is the vector of derivatives, ie. the RHS's of ODE's
# y[]:      dependent variables, initially = initial values 
# a:        start independent variable, end-point returned 
# h:        step-size to try, last used is returned 
# da:       total increment in independent variable 
# hmx:      maximum step allowed ever 
# iter:     maximum number of iterations allowed 
#
# rkf returns a success code as given below:
#      returned    status
#        +ve    success, and return value = iteration count
#        -1     exceeded iteration limit
#        -2     bad values supplied for atol, rtol, h or da.
#        -3     h, the step-size, has become too small for machine precision
# In addition, a, h, and y[...] are returned. a is the end-point, and would
# be the next starting-point to continue integration. h is the step-size last
# used, and is usually a valid guess of what step to use on next interval.
#
# The return of status, a and h is in a tuple: (status,a,h)

  
  hmx= da;
 
  if (rtol < 0) | (atol < 0) | (hmx <= 0) | \
                    ( (rtol == 0) & (atol == 0) ): return (-2,a,h)
  
  b = a + da                          # end point of the integration
  if abs(da) <= 10.0 * eps * max(abs(a), abs(b)) :
    return (-3,a,h)                   # da is too small.
  
  hmax = min(hmx, abs(da))
  if abs(h) <= 10.0 * eps * abs(a) : h = hmax
  kount = 0                           # zero function counter
  lasth = 0                           # not last step

  #-------------------------------------------------------------------------
  while 1 :                              # continue with new steps
    h = sign(min(abs(h), hmax), da)      # transfer sign of da to h
    if abs(b - a) <= 1.25 * abs(h) :     # close enough to leap
      h = b - a                          # size of last step required
      lasth = 1
    k1 = F(y, params)                          # call the RHS function
    kount = kount + 1
  
    #-----------------------------------------------------------------------
    while 1 :                            # no good, try a new step size
      tmp = y + a21 * h * k1
       
      arg = a + c2 * h 
      k2 = F(tmp, params)
      tmp = y + h * (a31 * k1 + a32 * k2)
      
      arg = a + c3 * h 
      k3 = F(tmp, params)      
      tmp = y + h * (a41 * k1 + a42 * k2 + a43 * k3)
      
      arg = a + c4 * h 
      k4 = F(tmp, params)      
      tmp = y + h * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
      
      arg = a + c5 * h 
      k5 = F(tmp, params)
      tmp = y + h * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
     
      arg = a +  c6 * h
      k6 = F(tmp, params)
      tmp = y + h * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5)
      maxerr = h * ( (bb1-b1)*k1 + (bb3-b3)*k3 + (bb4-b4)*k4 + (bb5-b5)*k5 + bb6*k6)
      
      maxerr = abs(maxerr) / (rtol * abs(tmp) + atol) 
      
      #y = tmp
      #return(1,tmp, h)
      if maxerr <= 1.0 :                  
        y = tmp        
        a = a + h                        
        
      
        if lasth == 1 :
          return (1,y,h, maxerr)             
        if maxerr > 6.5536e-4 :
          h = h * (0.8/sqrt(sqrt(maxerr)))      
        else : h = h * 5.0 
      
        return (1,y,h,maxerr) 
        break                                   
  
      if maxerr < 4096 :            
        h = h * (0.8/sqrt(sqrt(maxerr)))         # errors are 5th order
      else : h = h * 0.1 
      
      if abs(h) <= 10 * eps * max(abs(a), abs(b)) :
        return (-3,y,h, maxerr)                         # h has become too small.
      
      lasth = 0                       # step too large, is no longer last step
    #-----------------------------------------------------------------------
  #-------------------------------------------------------------------------  
  
  
  
  
  
  
  

#print('%0.12e'%eps)

#Read the (p,q,r) data file the 1st line contains the values of the parameters : 
# alpha, m, n, lambda, a, b 
# first read the parameters
#params=np.zeros(8)
#fp = open('myshbndPQR.dat')
#line = fp.readline(); line = fp.readline();  prm=line.split();
#print(line)
#alf=params[3]=float(prm[0])
#m=params[4]=float(prm[1])
#n=params[5]=float(prm[2])
#lam=params[6]=float(prm[3])
#aa=params[7]=float(prm[4])
#bb=float(prm[5])

params=np.zeros(8)
mypar = np.load('my3par.npy') 
alf=params[3]=float(mypar[0])
m=params[4]=float(mypar[1])
n=params[5]=float(mypar[2])
lam=params[6]=float(mypar[3])

myeps0 = mypar[6]

X01 = np.zeros(4)
X02 = np.zeros(4)
X03 = np.zeros(4)
X04 = np.zeros(4)
X11 = np.zeros(4)
X12 = np.zeros(4)
X14 = np.zeros(4)
X13 = np.zeros(4)

X01 = mypar[21:25]
X02 = mypar[25:29]
X03 = mypar[29:33]
X11 = mypar[33:37]
X12 = mypar[37:41]
X14 = mypar[41:45]

X13 = mypar[45:49]
X04 = mypar[49:53]

D = 1.0+2*alf-m-n
a0,a1 = (2+2*alf-n)   / D   , (2+2*alf)       / D 
b0,b1 = (1+m)       / D   , (1+m+n)       / D 

aa = params[7]=a0 + a1*lam
bb = b0 + b1*lam

#print(alpha,m,n,lam,aa,bb)
# now read the t,p,q,r data
#TPQR=np.loadtxt('myshbndPQR.dat', skiprows=2);
TPQR=np.loadtxt('myshbndPQR.dat');
N=len(TPQR[:,0]);
T=np.zeros(N); P=np.zeros(N); Q=np.zeros(N); R=np.zeros(N); S=np.zeros(N); 
T=TPQR[:,0]
P=TPQR[:,1]
Q=TPQR[:,2]
R=TPQR[:,3]


t= T[0]; Tfinal = T[N-1]
nstep = 0

# initial value for s 
DD=1.0+2.0*alf-m-n
r0=(2.0+2.0*alf-n)/DD + (2.0+2.0*alf)*lam/DD
S[0]=(1.0+m+n)/(1.0+alf) - n/((1.0+alf)*r0);

#myc01= 9.8572084967e-01
#myc02= 1.6838766735e-01
#myc03= 1.5512729904e-06 
#myeps0 = 1.0488740383e-08

S[0]=(1.0+m+n)/(1.0+alf) - n/((1.0+alf)*r0) + myeps0*mypar[8]*X01[3] + myeps0*mypar[9]*X02[3] + myeps0*mypar[10]*X03[3]

while (t < Tfinal ):
    dt = T[nstep+1] - T[nstep]
    ts = T[nstep]
    te = T[nstep+1]
    params[0] = P[nstep]
    params[1] = Q[nstep]
    params[2] = R[nstep]

    ss=S[nstep]; dtc = dt
    k, tmp, dtc, mxerr = rkf45(F,ss,ts,dtc,dt, params)   
    if (abs(dtc-dt)>10*eps): print("Warning : Next time has not be reached!!!")
        
    print(t, nstep,ss, tmp, mxerr)
    nstep = nstep + 1
    t = T[nstep]
    S[nstep] = tmp

#myt = T/20.0 + 0.5  

MAT = np.zeros((4,4))
MAT = np.transpose([X11,X12,X13,X14])
myf = np.transpose([P[-1],Q[-1]-1.0,R[-1]-0.7619047619047623,S[-1]-0.46875])/mypar[7]
myc = np.linalg.solve(MAT,myf)
    

  
np.savetxt('myshbnd.dat', np.transpose([T, P, Q , R, S]), fmt='%23.15E')
    
    
#fig, ax = plt.subplots()
#
#ax.plot(T,S)   
    
    
    
    
    
    

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(P, Q, R, marker='o')
#ax.set_xlabel('P')
#ax.set_ylabel('Q')
#ax.set_zlabel('R')
#
#
#fig = plt.figure()
#ax1 = fig.add_subplot(111, projection='3d')
#ax1.scatter(P, Q, S, marker='o')
#ax1.set_xlabel('P')
#ax1.set_ylabel('Q')
#ax1.set_zlabel('S')

