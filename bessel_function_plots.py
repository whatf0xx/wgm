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
n0 = 0 #Lower bound to generate
n = 4 #Upper bound to generate ('Hankel crashes at n=5?')

J = [spl.jv(n0, x)] #Bessel function, order 0; argument x
for i in range(n0, n):
    J.append(spl.jv(i+1, x)) #Create Bessel functions, order 1-->n; argument x

H = [spl.hankel1(n0, x)] #Hankel function, order 0; argument x
for j in range(n0, n):
    H.append(spl.hankel1(j+1, x)) #Create Bessel functions, order 1-->n; argument x

#%% plots
plt.subplot(2, 1, 1) #Bessel functions
plt.title("Bessel functions")
plt.grid(True)
plt.vlines(0, -0.5, 1.1)
plt.hlines(0, -1, 10)
for k in range(n0, n+1):
    plt.plot(x, J[k-n0], label="n=%d"% k)
plt.legend()
    
plt.subplot(2, 1, 2) #Hankel functions
plt.title("Hankel functions of the first kind")
plt.legend()
plt.grid(True)
plt.vlines(0, -0.5, 1.1)
plt.hlines(0, -1, 10)
for l in range(n0, n+1):
    plt.plot(x, H[l-n0], label="n=%d"% l)
plt.legend()
