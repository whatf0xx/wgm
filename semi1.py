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

plt.subplot(311)
plt.title("Cavity measurements of sound decay, $1/r$ fit")
n = 10001
xps = sp.linspace(0.36, 1.1, n)
def amp(x, a, b, c):
    return a/(x-b) + c
fit = curve_fit(amp, run1, amp1, p0=[1,0,0])
data_fit = amp(xps, *fit[0])
snd = sp.sqrt(sp.diag(fit[1]))
perc = abs(snd/fit[0]*100)
plt.plot(xps, data_fit, color='red', label=u"a=%.3f, b=%.3f, c=%.3f; $\u03C3_a$=%.3f, $\u03C3_b$=%.3f, $\u03C3_c$=%.3f"% (*fit[0], *snd))
plt.plot(run1, amp1, 'bx', label="percentage errors: a=%.1f, b=%.1f, c=%.1f"% (perc[0], perc[1], perc[2]))
plt.xlabel("Distance/m")
plt.ylabel("Amplitude/mV")
plt.legend()

plt.subplot(312)
plt.title("Cavity measurements of sound decay $1/r^2$ fit")
n = 10001
xps = sp.linspace(0.36, 1.1, n)
def amp(x, a, b, c):
    return a/(x-b)**2 + c
fit = curve_fit(amp, run1, amp1, p0=[1,0,0])
data_fit = amp(xps, *fit[0])
snd = sp.sqrt(sp.diag(fit[1]))
perc = abs(snd/fit[0]*100)
plt.plot(xps, data_fit, color='red', label=u"a=%.3f, b=%.3f, c=%.3f; $\u03C3_a$=%.3f, $\u03C3_b$=%.3f, $\u03C3_c$=%.3f"% (*fit[0], *snd))
plt.plot(run1, amp1, 'bx', label="percentage errors: a=%.1f, b=%.1f, c=%.1f"% (perc[0], perc[1], perc[2]))
plt.xlabel("Distance/m")
plt.ylabel("Amplitude/mV")
plt.legend()
plt.tight_layout()

plt.subplot(313)
plt.title("Cavity measurements of sound decay, $1/r + 1/r^2$ fit")
n = 10001
xps = sp.linspace(0.36, 1.1, n)
def amp(x, a, b, c, d, e):
    return a/(x-b) + c/(x-d)**2 + e
fit = curve_fit(amp, run1, amp1, p0=[1,1,1,1,1])
data_fit = amp(xps, *fit[0])
snd = sp.sqrt(sp.diag(fit[1]))
perc = abs(snd/fit[0]*100)
plt.ylim(0, 1.5)
plt.plot(xps, data_fit, color='red', label=u"a=%.3f, b=%.3f, c=%.3f, d=%.3f, e=%.3f; $\u03C3_a$=%.3f, $\u03C3_b$=%.3f, $\u03C3_c$=%.3f, $\u03C3_d$=%.3f, $\u03C3_e$=%.3f"% (*fit[0], *snd))
plt.plot(run1, amp1, 'bx', label="percentage errors: a=%.1f, b=%.1f, c=%.1f, d=%.1f, e=%.1f"% (perc[0], perc[1], perc[2], perc[3], perc[4]))
plt.xlabel("Distance/m")
plt.ylabel("Amplitude/mV")
plt.legend()