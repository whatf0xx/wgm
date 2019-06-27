# -*- coding: utf-8 -*-
"""
Created on Sun Jun 23 18:46:58 2019

@author: whatf
"""

import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

plt.subplot(2,1,1)
plt.title("Initial free space measurements with 'bump'")
d2, a2 = sp.loadtxt("freespace 2.txt", skiprows=1, unpack=True)
dist2 = d2/100
def amp(x, a, c, d):
    return a/(x-c)**2 + d
x = sp.linspace(0.04, 0.55, 10001)
opt, cov = curve_fit(amp, dist2, a2, [0.5, -0.25, -4])
snd = sp.sqrt(sp.diag(cov))
plt.xlabel('distance/m')
plt.ylabel('amplitude/V')
perc = snd/opt*100
plt.errorbar(dist2, a2, fmt='g.', label="percentage errors: a=%.1f, b=%.1f, c=%.1f"% (perc[0], perc[1], perc[2]))
plt.plot(x, amp(x, *opt), color='pink', label=u"a=%.3f, b=%.3f, c=%.3f; $\u03C3_a$=%.3f, $\u03C3_b$=%.3f, $\u03C3_c$=%.3f"% (*opt, *snd))
plt.legend()

plt.subplot(2,1,2)
plt.title("Later free space measurements elevated and with anti-echo foam")
da, am = sp.loadtxt("freegood.txt", unpack=True)
dista = da/100
def amp(x, a, c, d):
    return a/(x-c)**2 + d
x = sp.linspace(0.06, 1, 10001)
opt, cov = curve_fit(amp, dista, am, [0.5, -0.25, -4])
snd = sp.sqrt(sp.diag(cov))
plt.xlabel('distance/m')
plt.ylabel('amplitude/V')
perc = snd/opt*100
plt.errorbar(dista, am, fmt='b.', label="percentage errors: a=%.1f, b=%.1f, c=%.1f"% (perc[0], perc[1], perc[2]))
plt.plot(x, amp(x, *opt), color='red', label=u"a=%.3f, b=%.3f, c=%.3f; $\u03C3_a$=%.3f, $\u03C3_b$=%.3f, $\u03C3_c$=%.3f"% (*opt, *snd))
plt.legend()

plt.tight_layout()