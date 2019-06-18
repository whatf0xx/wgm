# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 11:41:30 2019

@author: hf4218
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
ks = [0.768, 1.03, 1.28, 1.52, 1.76, 1.99, 2.22, 2.45]

    #for this simulation, k is given and based on that several Bessel functions are calculated and superimposed.
fstr = input("input arbitrary frequency in Hz: ")
f = int(fstr)
k = f*2*sp.pi/c_air #calculate k based off input frequency

    #starting at m=1, work up to the highest order mode that fits
limit_reached = False
m = 7
wgm = r
while(limit_reached == False):
    
    if(k < ks[m]):
        limit_reached = True
    
    int_wgm_rho = spl.jv(m, n*r*k)
    ext_wgm_rho = spl.hankel1(m, r*k)
    B = spl.jv(m, n*a*k)/spl.hankel1(m, a*k)
    wgm_phi = sp.cos(m*phi) + sp.sin(m*phi)*sp.sqrt(-1)

    heavymap = sp.heaviside((r/a)-1, 1)
    if(m == 7):
        wgm = sp.real(wgm_phi*int_wgm_rho*(1-heavymap) + B*wgm_phi*ext_wgm_rho*heavymap)
    else:
        wgm = wgm + sp.real(wgm_phi*int_wgm_rho*(1-heavymap) + B*wgm_phi*ext_wgm_rho*heavymap)
        
    m += 1

#%%draw on circle
heavyp = sp.heaviside(r-0.995*a, 1)
heavyn = sp.heaviside(1.005*a-r, 1)
xlab = sp.round_(x, 2)
ylab = sp.round_(y, 2)
circle = pd.DataFrame(heavyp*heavyn, index=ylab, columns=xlab)

#%% make heatmap of function
g=sns.heatmap(wgm, cmap="RdBu_r", xticklabels=True, yticklabels=True, square=True)
circle_palette = [(0xFF/0xFF, 0xFF/0xFF, 0xFF/0xFF, 0.001), (0xD1/0xFF, 0xEC/0xFF, 0x9C/0xFF, 1)]
cmap = mpl.colors.ListedColormap(circle_palette)
c=sns.heatmap(circle, cmap=cmap, xticklabels=100, yticklabels=100, cbar=False)
g.set_title("f = %f Hz"% f)
g.set_xlabel("X (m)")
g.set_ylabel("Y (m)")
