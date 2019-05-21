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
x = sp.linspace(-1.5, 1.5, 1000)
y = sp.linspace(-1.5, 1.5, 1000)
X, Y = sp.meshgrid(x, y)

    #convert to plane polars
r = sp.sqrt(X**2 + Y**2)
phi = sp.arctan2(Y, X)

    #set constants
c_air = 343 #speed of sound in air (m/s)
c_mat = 3500 #speed of sound in barrier material/outside cavity (m/s). Here the value used is for concrete
n = c_mat/c_air #cavity refractive index
f = 700 #frequency of sound waves (Hz)
omega = f*2*sp.pi
k = omega/c_mat #wavenumber (/m)
a = 1 #arbitrary cavity radius (m)
m = 9

    #evaluate functions
int_wgm_rho = spl.jv(m, n*r*k)
ext_wgm_rho = spl.hankel1(m, r*k)
B = spl.jv(m, n*a*k)/spl.hankel1(m, a*k)
wgm_phi = sp.cos(m*phi) + sp.sin(m*phi)*sp.sqrt(-1)
heavymap = sp.heaviside((r/a)-1, 1)

    #draw on circle
x_circ = sp.linspace(-a, a, 1000)
y_circ_pos = sp.sqrt(a**2-x_circ**2)
y_circ_neg = -sp.sqrt(a**2-x_circ**2)

#%% calculate values
wgm = sp.real(wgm_phi*int_wgm_rho*(1-heavymap) + B*wgm_phi*ext_wgm_rho*heavymap)

#%% make heatmap of function
sns.heatmap(wgm)
plt.plot(x_circ, y_circ_pos, color='green')
plt.plot(x_circ, y_circ_neg, color='green')
plt.show()
