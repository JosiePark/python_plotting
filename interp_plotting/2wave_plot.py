#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import sys

sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import square_grid_plot

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/INTERP/TEST_CASES/TWO_WAVE/'
fig_dir = '/home/Project/PhD/PYTHON_FIGURES/INTERP/2WAVE/'
ii = 512

# READ PSI, U AND V #

psi_file = open(home_dir + 'NO_RK4/2wave.dat','r')
psi = [ord[i] for i in psi_file.read()]
print(psi)
print(psi[0])
print(float(psi[-1]))
sys.exit()

u_file = home_dir + 'RK4/2wave_u.dat'
u = np.fromfile(u_file,dtype=float)

v_file = home_dir + 'RK4/2wave_v.dat'
v = np.fromfile(v_file,dtype=float)

# PLOT PSI, U AND V #

nrows = 1
ncols = 3
top_space = .05
bottom_space = 0.08
right_space = 0.05
left_space = 0.05
ver_space = 0.
hor_space = 0.05
fig_width = 8.27

ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

xmin,xmax = 0,520
xx,yy = np.mgrid[xmin:xmax:512j,xmin:xmax:512j]

k = 0

for i in range(ncols):
	ax[k].contourf(xx,yy,psi,cmap = cm.jet)
	
plt.show()

