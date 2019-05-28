# -*- coding: utf-8 -*-
"""
Created on Tue May 28 10:28:08 2019

@author: @team_effort
"""

#%% take 1, a bit of a mare
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

distance, amplitude = sp.loadtxt("freespace_sound.txt", skiprows=1, unpack=True)
d = distance/100 #cm to m
plt.plot(d, amplitude, 'bx') #initial readings were so bad lol
I = amplitude*d**2
plt.plot(d, I)

#%% take 2
d2, a2 = sp.loadtxt("freespace 2.txt", skiprows=1, unpack=True)
dist2 = d2/100
plt.plot(dist2, a2)
I2 = a2*dist2**2
plt.plot(dist2, I2)

#%% take 3 - actually quite good
d3, a3 = sp.loadtxt("free3.txt", unpack=True)
dist3 = d3/100
plt.plot(dist3, a3, 'bx')
I3 = 1/sp.sqrt(a3)
#plt.plot(dist3, I3)
def amp(x, a, c, d):
    return a/(x-c)**2 + d
x = sp.linspace(-0.1, 0.5, 500)
y = amp(x, 0.5, -0.25, -3)
opt, cov = curve_fit(amp, x, y, [0.5, -0.25, -4])
plt.plot(x, amp(x, *opt))
plt.plot(x, y) 
#grad, inte = sp.polyfit(dist3, I3, 1)
#fit = (dist3*grad)**-2
#plt.plot(dist3, fit)