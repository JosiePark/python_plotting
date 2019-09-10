#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import sys
import scipy.io

sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import square_grid_plot

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/INTERP/TEST_CASES/TWO_WAVE/'
fig_dir = '/home/Project/PhD/PYTHON_FIGURES/INTERP/2WAVE/'
ii = 512

# READ PSI, U AND V #

psi_file = home_dir + 'NO_RK4/2wave_psi.nc'
psi_data = Dataset(psi_file,'r')
psi = psi_data.variables['psi'][:]

u_file = home_dir + 'NO_RK4/2wave_u.nc'
u_data = Dataset(u_file,'r')
u = u_data.variables['u'][:]

v_file = home_dir + 'NO_RK4/2wave_v.nc'
v_data = Dataset(v_file,'r')
v = v_data.variables['v'][:]

# PLOT PSI, U AND V #

nrows = 1
ncols = 3
top_space = .15
bottom_space = 0.2
right_space = 0.02
left_space = 0.08
ver_space = 0.
hor_space = 0.
fig_width = 8.27

ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

xmin,xmax = 0,520
xx,yy = np.mgrid[xmin:xmax:512j,xmin:xmax:512j]
	
ax[0].contourf(xx,yy,psi,cmap = cm.jet)
ax[0].set_title('Psi')
ax[0].set_ylabel('Y (km)')
ax[0].set_xlabel('X (km)')
ax[1].contourf(xx,yy,u,cmap = cm.jet)
ax[1].set_title('U')
ax[1].set_xlabel('X (km)')
ax[1].set_yticklabels([])
ax[2].contourf(xx,yy,v,cmap=cm.jet)
ax[2].set_title('V')
ax[2].set_xlabel('X (km)')
ax[2].set_yticklabels([])
	
plt.savefig('2wave.png')

