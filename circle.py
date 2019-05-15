# -*- coding: utf-8 -*-
"""
Created on Mon May 13 15:12:09 2019

@author: hf4218
"""

#%% libraries
import scipy as sp
import matplotlib.pyplot as plt

#%% create and draw the circle
x_points = sp.linspace(-1, 1, 100000)
circle_positive = sp.sqrt(1-x_points**2)
circle_negative = (-1)*sp.sqrt(1-x_points**2)

plt.plot(x_points, circle_positive)
plt.plot(x_points, circle_negative)

#%% initialise a start point at the base of the circle
line_x = [[0]]
line_y = [[-1]]

#%% calculate line gradient for beam
n = 8
theta = [sp.pi/n]
m = [sp.tan(theta[0])]
c = [-1]
#%% begin loop to draw lines
p = 0
while(p<n):
    q=0
    while(sp.sqrt(line_x[p][q]**2+line_y[p][q]**2)<=1):
        if((theta[p]>sp.pi/2 and theta[p]<sp.pi*3/2) or (theta[p]<-sp.pi/2 and theta[p]>-sp.pi*3/2)):
            line_x[p].append(line_x[p][q]-0.00001)
        else:
            line_x[p].append(line_x[p][q]+0.00001)
        line_y[p].append(line_x[p][q+1]*m[p]+c[p])
        q = q+1
    del line_x[p][q]
    del line_y[p][q]
    q = q-1
    plt.plot(line_x[p], line_y[p])
    
    theta.append(sp.pi-theta[p]+2*sp.arctan2(line_y[p][q], line_x[p][q]))
    m.append(sp.tan(theta[p+1]))
    c.append(line_y[p][q]-m[p+1]*line_x[p][q])
    
    line_x.append([line_x[p][q]])
    line_y.append([line_y[p][q]])
    
    p = p+1
