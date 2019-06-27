# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 13:59:03 2019

@author: whatf
"""

import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

d,a = sp.loadtxt("radial measurements.txt", skiprows=1, unpack=True)
d = d/100
plt.plot(d, a, 'bx')
plt.ylim(100, 250)

average = sp.sum(a)/len(a)
plt.hlines(average, 0, 0.4, color='red', label="arithmetic mean of points")

def fit(x, a, b, c):
    return a*(x-b)**2 +c

fitted = curve_fit(fit, d, a, p0=[0.1,0.2, 0])
snd = sp.sqrt(sp.diag(fitted[1]))
x = sp.linspace(0, 0.4, 101)
plt.plot(x, fit(x, *fitted[0]), color='green', label="parabolic fit")

plt.title("Variations in amplitude with respect to radial distance")
plt.xlabel("Distance from the cavity edge/m")
plt.ylabel("Amplitude/mV")
plt.legend()