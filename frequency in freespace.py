# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 13:21:11 2019

@author: whatf
"""

import scipy as sp
import matplotlib.pyplot as plt

    #what a horrible characteristic
frequency, amplitude = sp.loadtxt("frequency free space.txt", skiprows=1, unpack=True)
plt.plot(frequency, amplitude, 'bx')
plt.title("Amplitude against frequency characteristic for the setup used")
plt.xlabel("Frequency/Hz")
plt.ylabel("Amplitude/mV")

frequency, cavamp = sp.loadtxt("cavity frequency.txt", skiprows=1, unpack=True)
adjamp = cavamp - amplitude

#plt.plot(frequency, adjamp, 'ro')