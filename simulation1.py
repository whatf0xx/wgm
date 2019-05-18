# -*- coding: utf-8 -*-
"""
Created on Sat May 18 18:23:48 2019

@author: whatf
"""

#%% libraries
import scipy as sp
import scipy.special as spl
import matplotlib.pyplot as plt
import seaborn as sns

#%% define variables and functions
    #grid
x = sp.linspace(-2, 2, 1000)
y = sp.linspace(-2, 2, 1000)
X, Y = sp.meshgrid(x, y)

    #convert to plane polars
r = sp.sqrt(X**2 + Y**2)
phi = sp.tan(Y/X)

    #set constants
c_air = 343 #speed of sound in air (m/s)
c_mat = 3500 #speed of sound in barrier material/outside cavity (m/s). Here the value used is for concrete
n = c_mat/c_air #cavity refractive index
f = 2000 #frequency of sound waves (Hz)
omega = f*2*sp.pi
k = omega/c_mat #wavenumber (/m)
a = 1 #arbitrary cavity radius (m)
m = 5
    #evaluate functions
int_wgm_rho = spl.jv(m, n*r*k)
ext_wgm_rho = spl.hankel1(m, r*k)
B = spl.jv(m, n*a*k)/spl.hankel1(m, a*k)
wgm_phi = sp.cos(m*phi)

#%% calculate values
for i in range (0,len(X)):
    for j in range (0, len(Y)):
        if(r[i][j]<a):
            wgm = wgm_phi*int_wgm_rho
        else:
            wgm = B*wgm_phi*ext_wgm_rho

#%% make heatmap of function
ax = sns.heatmap(wgm, linewidth=0.5)
plt.show()