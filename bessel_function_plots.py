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
x = sp.linspace(1, 10, 5000) #Note for higher order Hankel functions exclude x values close to the origin
n0 = 48 #Lower bound to generate
n = 49 #Upper bound to generate
a = 10.2 #Arbitrary constant corresponding to refractive index

J = [spl.jv(n0, x)] #Bessel function, order n0; argument x
H = [spl.hankel1(n0, x)] #Hankel function, order n0; argument x
Ja = [spl.jv(n0, x*a)] #Bessel function, order n0; argument a*x
Jdash = [spl.jvp(n0, x, n=1)] #Bessel function first derivative, order n0; argument x
Hdash = [spl.h1vp(n0, x, n=1)] #Hankel function first derivative, order n0; argument x
Jadash = [spl.jvp(n0, x*a, n=1)] #Bessel function first derivative, order 0; argument a*x
for i in range(n0, n):
    J.append(spl.jv(i+1, x)) #Create Bessel functions, order n0+1-->n; argument x
    H.append(spl.hankel1(i+1, x)) #Create Hankel functions, order n0+1-->n; argument x
    Ja.append(spl.jv(i+1, x*a)) #Create Bessel functions, order n0+1-->n; argument a*x
    Jdash.append(spl.jvp(i+1, x, n=1)) #Create Bessel first derivatives, order n0+1-->n; argument x
    Hdash.append(spl.h1vp(i+1, x, n=1)) #Create Hankel first derivatives, order n0+1-->n; argument x
    Jadash.append(spl.jvp(i+1, x*a, n=1)) #Create Bessel first derivatives, order 1-->n; argument a*x

eigroot = [Jadash[0]/Ja[0] - Hdash[0]/(a*H[0])]
for i in range(0, n-n0):
    eigroot.append(Jadash[i+1]/Ja[i+1] - Hdash[i+1]/(a*H[i+1]))
#%% plots
plt.subplot(3, 1, 1) #Bessel functions
plt.title("Bessel functions")
plt.grid(True)
plt.vlines(0, -0.3, 0.3)
plt.hlines(0, -1, 15)
for i in range(n0, n+1):
    plt.plot(x, J[i-n0], label="n=%d"% i)
plt.legend()
    
plt.subplot(3, 1, 2) #Hankel functions
plt.title("Hankel functions of the first kind")
plt.legend()
plt.grid(True)
plt.vlines(0, -0.3, 0.3)
plt.hlines(0, -1, 15)
for i in range(n0, n+1):
    plt.plot(x, sp.real(H[i-n0]), label="n=%d"% i)
plt.legend()

plt.subplot(3, 1, 3) #Eigenfrequency solutions
plt.title("Newton-Raphson solutions for eigenmodes")
plt.legend()
plt.grid(True)
plt.vlines(0, -0.3, 0.3)
plt.hlines(0, -1, 15)
for i in range(n0, n+1):
    plt.plot(x, sp.real(eigroot[i-n0]), label="n=%d"% i)
plt.legend()

plt.tight_layout()