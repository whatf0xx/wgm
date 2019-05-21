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
x = sp.linspace(1, 15, 5000) #Note for higher order Hankel functions exclude x values close to the origin
n0 = 7 #Lower bound to generate
n = 9 #Upper bound to generate

J = [spl.jv(n0, x)] #Bessel function, order 0; argument x
for i in range(n0, n):
    J.append(spl.jv(i+1, x)) #Create Bessel functions, order 1-->n; argument x

H = [sp.real(spl.hankel1(n0, x))] #Hankel function, order 0; argument x
for j in range(n0, n):
    H.append(sp.real(spl.hankel1(j+1, x))) #Create Hankel functions, order 1-->n; argument x

#%% plots
plt.subplot(2, 1, 1) #Bessel functions
plt.title("Bessel functions")
plt.grid(True)
plt.vlines(0, -0.3, 0.3)
plt.hlines(0, -1, 15)
for k in range(n0, n+1):
    plt.plot(x, J[k-n0], label="n=%d"% k)
plt.legend()
    
plt.subplot(2, 1, 2) #Hankel functions
plt.title("Hankel functions of the first kind")
plt.legend()
plt.grid(True)
plt.vlines(0, -0.3, 0.3)
plt.hlines(0, -1, 15)
for l in range(n0, n+1):
    plt.plot(x, H[l-n0], label="n=%d"% l)
plt.legend()
