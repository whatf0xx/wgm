# -*- coding: utf-8 -*-
"""
Created on Mon May 20 12:00:18 2019

@author: whatf
"""

#%% import libraries
import scipy as sp
import scipy.special as spl
import matplotlib.pyplot as plt

#%% create axes
x = sp.linspace(-1, 10, 50000)
J = [spl.jv(0, x)] #Bessel function, order 0; argument x
n = 4
for i in range(0, n):
    J.append(spl.jv(i+1, x)) #Create Bessel functions, order 1-->n; argument x

H = [spl.hankel1(0, x)] #Hankel function, order 0; argument x
for j in range(0, n):
    H.append(spl.hankel1(j+1, x)) #Create Bessel functions, order 1-->n; argument x

#%% plots
plt.subplot(2, 1, 1) #Bessel functions
plt.title("Bessel functions")
plt.grid(True)
plt.vlines(0, -0.5, 1.1)
plt.hlines(0, -1, 10)
for k in range(0, n+1):
    plt.plot(x, J[k], label="n=%d"% k)
plt.legend()
    
plt.subplot(2, 1, 2) #Hankel functions
plt.title("Hankel functions of the first kind")
plt.legend()
plt.grid(True)
plt.vlines(0, -0.5, 1.1)
plt.hlines(0, -1, 10)
for l in range(0, n+1):
    plt.plot(x, H[l], label="n=%d"% l)
plt.legend()
