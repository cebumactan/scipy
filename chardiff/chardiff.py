# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 18:46:22 2016

@author: leem
"""
import numpy as np
import matplotlib.pyplot as plt

# import data
mydata = np.load('./mydata.npz')
Hplus = mydata['arr_0']
Voltage = mydata['arr_2']
t = mydata['arr_3']
pH = mydata['arr_4']

figure(1)
plt.plot(t,pH)

figure(2)
plt.plot(t,Hplus)
