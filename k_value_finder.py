# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 18:18:27 2019

@author: whatf
"""

#%% libraries
import scipy as sp
import scipy.special as spl
#import pandas as pd

#%% variables
c_air = 343 #speed of sound in air (m/s)
c_mat = 3000 #speed of sound in barrier material/outside cavity (m/s). Here the value used is for concrete
n = c_mat/c_air #cavity refractive index
a = 0.57 #arbitrary cavity radius (m)

while(True):
    m_str = input("Choose the m value of the Bessel function (or use '-1' to exit): ")
    if(m_str == '-1'):
        break
    k_str = input("Guess a corresponding k value: ")

    m = int(m_str)
    k = float(k_str)

    done = False
    counter = 0

#%% perform NR method
    while(done==False):
        J = spl.jv(m, n*a*k) #Bessel function at boundary (a)
        Jdash = spl.jvp(m, n*a*k, n=1) #first derivative
        Jdodash = spl.jvp(m, n*a*k, n=2) #second
        H = spl.hankel1(m, a*k) #Hankel function at boundary
        Hdash = spl.h1vp(m, a*k, n=1) #first derivative
        Hdodash = spl.h1vp(m, a*k, n=2) #second
        
        NRfunc = Jdash/J - Hdash/(n*H) #compute Newton-Raphson approximation for eigenfrequency
        NRfuncdash = n*a*(Jdodash*J - 2*Jdash)/(J**2) - a*(Hdodash*H - 2*Hdash)/(n*H**2)
        k_prev = k
        k = k - NRfunc/NRfuncdash
        counter += 1
        if(abs(k-k_prev) < 0.0001):
            done=True

#%% print results
    f = k.real*c_air/(2*sp.pi)
    print("f = %.1f Hz"% f)