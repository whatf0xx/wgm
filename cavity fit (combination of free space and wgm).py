import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

d3, a3 = sp.loadtxt("cavity measurements.txt", unpack=True)
dist3 = d3/100

plt.plot(dist3, a3, 'bx')
n = 10001
xps = sp.linspace(0.36, 1.1, n)
def amp(x, a, b, c, d, e):
    return a/(x-b)**2 + c/(x-d) + e

fit = curve_fit(amp, dist3, a3, p0=[1,0,1,0,0])

data_fit = amp(xps, *fit[0])
plt.xlabel('distance/m')
plt.ylabel('amplitude/V')
utol = 2.5
ltol = 0
data_fit[data_fit>utol] = sp.nan
data_fit[data_fit<ltol] = sp.nan
plt.plot(xps, data_fit)
plt.plot(dist3, a3, 'x')
plt.grid()
plt.show()

#%%
print("fit variables for cavity used: a=%.5f, b=%.5f, c=%.5f, d=%.5f, e=%.5f"% (fit[0][0], fit[0][1], fit[0][2], fit[0][3],fit[0][4]))
snd = sp.sqrt(sp.diag(fit[1]))
print("these has standard deviations of a=%.4f, b=%.4f, c=%.4f, d=%.4f, e=%.4f"% (snd[0], snd[1], snd[2], snd[3], snd[4]))
