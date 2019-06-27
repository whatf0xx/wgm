# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 11:46:28 2019

@author: whatf
"""

import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

d3, a3 = sp.loadtxt("cavity measurements.txt", unpack=True)
dist3 = d3/100

run1 = dist3[:6]
amp1 = a3[:6]

run2 = dist3[6:]
amp2 = a3[6:]


plt.subplot(211)
plt.plot(run1, amp1, 'bx')
n = 10001
xps = sp.linspace(0.36, 1.1, n)
def amp(x, c, d, e):
    return c/(x-d) + e
fit = curve_fit(amp, run1, amp1, p0=[1,0,0])
data_fit = amp(xps, *fit[0])
plt.plot(xps, data_fit, color='red')

plt.subplot(212)
plt.plot(run2, amp2, 'gx')
fit = curve_fit(amp, run2, amp2, p0=[1,0,0])
data_fit = amp(xps, *fit[0])
plt.plot(xps, data_fit, color='pink')

#%%
print("fit variables for cavity used: a=%.5f, b=%.5f, c=%.5f, d=%.5f, e=%.5f"% (fit[0][0], fit[0][1], fit[0][2], fit[0][3],fit[0][4]))
snd = sp.sqrt(sp.diag(fit[1]))
print("these has standard deviations of a=%.4f, b=%.4f, c=%.4f, d=%.4f, e=%.4f"% (snd[0], snd[1], snd[2], snd[3], snd[4]))
plt.xlabel('distance/m')
plt.ylabel('amplitude/V')