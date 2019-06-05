# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 02:14:52 2019

@author: whatf
"""

import scipy as sp
import seaborn as sns

x = sp.linspace(-1.5, 1.5, 1000)
y = sp.linspace(-1.5, 1.5, 1000)
X, Y = sp.meshgrid(x, y)

r = sp.sqrt(X**2 + Y**2)
heavyp = sp.heaviside(r-0.995, 10)
heavyn = sp.heaviside(1.005-r, 10)
points = heavyp*heavyn
sns.heatmap(points) 