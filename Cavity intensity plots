import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

d3, a3 = sp.loadtxt("cavity measurements.txt", unpack=True)
dist3 = d3/100

plt.plot(dist3, a3, 'bx')
xps = sp.linspace(0.36, 1.1, 101)
def amp(x, a, b, c):
    return a/(x-b) + c

fit = curve_fit(amp, dist3, a3, p0=[-1,-1,1])

data_fit = amp(xps, *fit[0])
plt.xlabel('distance/m')
plt.ylabel('amplitude/V')
plt.plot(xps, data_fit)
plt.plot(dist3, a3, 'x')
plt.show()

print("fit variables used: a=%.3f, b=%.3f, c=%.3f."% (fit[0][0], fit[0][1], fit[0][2]))
snd = sp.sqrt(sp.diag(fit[1]))
print("these has standard deviations of a=%.3f, b=%.3f, c=%.3f."% (snd[0], snd[1], snd[2]))
#%%
theta = sp.arcsin(dist3/1.14)
plt.plot(theta, a3, 'bx')
xps = sp.linspace(0.32, 1.4, 101)
def amp(x, a, b, c):
    return a/(x-b) + c

fit = curve_fit(amp, theta, a3, p0=[-1,-1,1])

data_fit = amp(xps, *fit[0])
plt.xlabel('distance/m')
plt.ylabel('amplitude/V')
plt.plot(xps, data_fit)
plt.plot(theta, a3, 'x')
plt.show()

print("fit variables used: a=%.3f, b=%.3f, c=%.3f."% (fit[0][0], fit[0][1], fit[0][2]))
snd = sp.sqrt(sp.diag(fit[1]))
print("these has standard deviations of a=%.3f, b=%.3f, c=%.3f."% (snd[0], snd[1], snd[2]))
