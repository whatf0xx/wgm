# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:31:51 2019
@author: nellp
"""
#%%
import scipy as sp
import matplotlib.pyplot as plt
x_points = sp.linspace(-1, 1, 100000)
circle_positive = sp.sqrt(1-x_points**2)
circle_negative = (-1)*sp.sqrt(1-x_points**2)

plt.plot(x_points, circle_positive)
plt.plot(x_points, circle_negative)

line_x=[0]
line_y=[-1]
m=0.3
n=0

while(sp.sqrt(line_x[n]**2+line_y[n]**2)<=1):
    line_x.append(line_x[n]+0.0001)
    line_y.append(line_x[n+1]*m-1)
    n=n+1
plt.plot(line_x, line_y)

theta=sp.pi-sp.arctan(m)+2*sp.arctan2(line_y[n-1], line_x[n-1])
grad = sp.tan(theta)
inte = line_y[n-1]-grad*line_x[n-1]
line_x1=[line_x[n-1]]
line_y1=[line_y[n-1]]

x_p = sp.linspace(0.2, 0.8, 100)
y_p = x_p*line_y[n]/line_x[n]

n=0
while(sp.sqrt(line_x1[n]**2+line_y1[n]**2)<=1):
    line_x1.append(line_x1[n]+0.0001)
    line_y1.append(line_x1[n+1]*grad+inte)
    n=n+1
plt.plot(line_x1, line_y1)
plt.plot(x_p, y_p)
