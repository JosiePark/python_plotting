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

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/INTERP/TEST_CASES/STATIONARY_JET/'
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/INTERP/JET/'
ii = 512

# READ PSI, U AND V #

psi_file = home_dir + 'NO_RK4/QG.nc'
psi_data = Dataset(psi_file,'r')
psi = np.transpose(psi_data.variables['Stream Function'][0,0,:,:])


# PLOT PSI, U AND V #

nrows = 1
ncols = 1
top_space = .08
bottom_space = 0.1
right_space = 0.02
left_space = 0.1
ver_space = 0.
hor_space = 0.
fig_width = 8.27

ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

xmin,xmax = 0,520
xx,yy = np.mgrid[xmin:xmax:512j,xmin:xmax:512j]
	
c = ax[0].pcolor(xx,yy,psi[:,:],cmap = cm.jet)
ax[0].set_title('Psi (cm$^2$s$^{-1}$)')
ax[0].set_ylabel('Y (km)')
ax[0].set_xlabel('X (km)')
plt.colorbar(c,ax=ax[0])
	
plt.savefig(fig_dir + 'jet.png')

