import numpy as np
import matplotlib.pyplot as plt

u = np.loadtxt ('r_m')

x_arr = np.linspace (0, 2, u.shape[1])
y_arr = np.linspace (0, 2, u.shape[0])
x_arr, y_arr = np.meshgrid(x_arr, y_arr)


fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot_surface (x_arr, y_arr, u, cmap='viridis', edgecolor='none')
plt.show()