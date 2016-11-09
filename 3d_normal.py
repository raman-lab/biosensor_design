#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import bivariate_normal
from mpl_toolkits.mplot3d import Axes3D

#Parameters to set
mu_x = 6
sigma_x = np.sqrt(5)

mu_y = 7
sigma_y = np.sqrt(5)

#Create grid and multivariate normal
x = np.linspace(-10,10,500)
y = np.linspace(-10,10,500)
X, Y = np.meshgrid(x,y)
Z = bivariate_normal(X,Y,sigma_x,sigma_y,mu_x,mu_y)
Z = np.sin(x)*np.sin(y)
# Z = np.random.uniform(0,1,500)
Z = np.linspace(1,1,500)

#Make a 3D plot
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, Z,cmap='viridis',linewidth=0)
ax.set_axis_off()
ax.set_zlim(0,10)
plt.savefig('3d_normal.png', format='png', dpi=1000)