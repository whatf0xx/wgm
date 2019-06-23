# -*- coding: utf-8 -*-
"""
Created on Sat May 18 18:23:48 2019
@author: whatf
"""

#%% libraries
import scipy as sp
import scipy.special as spl
import matplotlib as mpl
import seaborn as sns
import pandas as pd

#%% define variables and functions
    #grid
x = sp.linspace(-1, 1, 1000)
y = sp.linspace(-1, 1, 1000)
X, Y = sp.meshgrid(x, y)

    #convert to plane polars
r = sp.sqrt(X**2 + Y**2)
phi = sp.arctan2(Y, X)

    #set constants
c_air = 343 #speed of sound in air (m/s)
c_mat = 3000 #speed of sound in barrier material/outside cavity (m/s). Here the value used is for concrete
n = c_mat/c_air #cavity refractive index
a = 0.57 #arbitrary cavity radius (m)
m = 9
k_trial = 2.7 #wavenumber (/m)
k_guess = k_trial #so the original value can be compared to the computed one
done = False
counter = 0
    #calculate eigenfrequencies
while(done==False):
    J = spl.jv(m, n*a*k_trial) #Bessel function at boundary (a)
    Jdash = spl.jvp(m, n*a*k_trial, n=1) #first derivative
    Jdodash = spl.jvp(m, n*a*k_trial, n=2) #second
    H = spl.hankel1(m, a*k_trial) #Hankel function at boundary
    Hdash = spl.h1vp(m, a*k_trial, n=1) #first derivative
    Hdodash = spl.h1vp(m, a*k_trial, n=2) #second
    
    NRfunc = Jdash/J - Hdash/(n*H) #compute Newton-Raphson approximation for eigenfrequency
    NRfuncdash = n*a*(Jdodash*J - 2*Jdash)/(J**2) - a*(Hdodash*H - 2*Hdash)/(n*H**2)
    k_prev = k_trial
    k_trial = k_trial - NRfunc/NRfuncdash
    counter += 1
    if(abs(k_trial-k_prev) < 0.0001):
        done=True

k = k_trial #the computed eigenfrequency, corresponding to the first root
f = k.real*c_air/2/sp.pi

    #evaluate functions
int_wgm_rho = spl.jv(m, n*r*k)
ext_wgm_rho = spl.hankel1(m, r*k)
B = spl.jv(m, n*a*k)/spl.hankel1(m, a*k)
wgm_phi = sp.cos(m*phi) + sp.sin(m*phi)*sp.sqrt(-1)
heavymap = sp.heaviside((r/a)-1, 1)

    #draw on circle
heavyp = sp.heaviside(r-0.995*a, 1)
heavyn = sp.heaviside(1.005*a-r, 1)
xlab = sp.round_(x, 2)
ylab = sp.round_(y, 2)
circle = pd.DataFrame(heavyp*heavyn, index=ylab, columns=xlab)

#%% calculate values
wgm = sp.real(wgm_phi*int_wgm_rho*(1-heavymap) + B*wgm_phi*ext_wgm_rho*heavymap)

#%% make heatmap of function
g=sns.heatmap(wgm, cmap="RdBu_r", square=True)
circle_palette = [(0xFF/0xFF, 0xFF/0xFF, 0xFF/0xFF, 0.001), (0xD1/0xFF, 0xEC/0xFF, 0x9C/0xFF, 1)]
cmap = mpl.colors.ListedColormap(circle_palette)
c=sns.heatmap(circle, cmap=cmap, xticklabels=100, yticklabels=100, cbar=False)
g.set_title(r"m = %d; k = %.3f+%.3fi m$^{-1}$; f = %.1f Hz"% (m, k.real, k.imag, f))
g.set_xlabel("X (m)")
g.set_ylabel("Y (m)")
