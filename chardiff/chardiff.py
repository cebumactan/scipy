# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 18:46:22 2016

@author: leem
"""
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

# import data
mydata = np.load('./mydata.npz')
Hplus = mydata['arr_0']
Voltage = mydata['arr_2']
t = mydata['arr_3']
pH = mydata['arr_4']
length = t.size
i_inj = 32540


# dimensional quantities
L = 4           #[mm]
T = 60          #[sec]
M = 4**3*1e-6   #[mol]
C = 1           #[mol/Liter] = 1000[mol/m^3] 
x0 = 38        #[mm], measuring point from the origin

D = np.zeros((length,))

# non-dimensional quantities
tt = (t-t[i_inj])/T
xx0 = x0/L
DD = np.zeros((length,))

# calculate the non-dimensonal diffusion constant at each time
DD[0]=0
DD[1]=0
for i in range(1,length-1):
    eta1 = np.log10( Hplus[i+1]*np.sqrt(tt[i+1]) )
    eta0 = np.log10( Hplus[i]  *np.sqrt(tt[i]  ) )
    DD[i+1] = ( ( tt[i+1]-tt[i] ) * ( xx0**2/(4*(tt[i]**2)) ) / np.log(10) ) / (eta1-eta0)

# dimensional diffusion constant
D = DD * L**2/T *1e-6

# Plot
fig=plt.figure(3,figsize=(45,8));plt.xlabel('time (min)'); plt.ylabel('diffusion constant(m^2 sec^(-1))');
plt.ticklabel_format(style='sci', axis='y', scilimits=(-9,-5))
plt.ylim([-1e-5,1e-5])
plt.plot(tt[50000:200000],D[50000:200000],'.')

plt.figure(1);plt.xlabel('time (min)'); plt.ylabel('pH'); plt.plot(np.log10(tt[i_inj:]),pH[i_inj:])
plt.figure(2);plt.xlabel('time (min)'); plt.ylabel('[H+] mol/Liter');plt.plot(tt[i_inj:],Hplus[i_inj:])

##a = np.zeros((length,))
#myD = 1e-2
#a = exp(-2**2/(4*myD*tt[1:]))/ sqrt(np.pi*D*tt[1:])
#figure(4);
#plt.plot(tt[1:],a) 
