#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import sys

sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import square_grid_plot

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/INTERP/TEST_CASES/TWO_WAVE/NO_RK4/'
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/INTERP/2WAVE/'
ii = 512

runs = 4
dt = np.zeros((runs,1))
for i in range(runs):
	dt[i] = 8640./(2**i)

nrows = 4
ncols = 3
left_space = 0.08
right_space = 0.02
bottom_space = 0.05
top_space = 0.05
hor_space = 0
ver_space = 0.05
fig_width = 8.27

ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

for i in range(runs):
	file_name = home_dir + '2wave_bicubic_dt_%i.nc' % dt[i]
	file_data = Dataset(file_name,'r')
	tmp = np.transpose(file_data.variables['Trajectories'][:])*520./ii
	ax[i].plot(tmp[0,0,:],tmp[1,0,:],'k')
	ax[i].set_title('Bicubic,\n dt = %i secs' % dt[i])
	if i == 3:
		ax[i].set_xlabel('X(km)')
	else:
		ax[i].set_xticklabels([])
	ax[i].set_ylabel('Y(km)')
	
	
	file_name = home_dir + '2wave_2Dcubic_dt_%i.nc' % dt[i]
	file_data = Dataset(file_name,'r')
	tmp = np.transpose(file_data.variables['Trajectories'][:])*520./ii
	ax[4+i].plot(tmp[0,0,:],tmp[1,0,:],'k')
	ax[4+i].set_title('2D cubic,\n dt = %i secs' % dt[i])
	if i == 3:
		ax[i+4].set_xlabel('X(km)')
	else:
		ax[i+4].set_xticklabels([])
	ax[i+4].set_yticklabels([])
	
	file_name = home_dir + '2wave_exact_dt_%i.nc' % dt[i]
	file_data = Dataset(file_name,'r')
	tmp = np.transpose(file_data.variables['Trajectories'][:])*520./ii
	ax[8+i].plot(tmp[0,0,:],tmp[1,0,:],'k')
	ax[8+i].set_title('No Interpolation,\n dt = %i secs' % dt[i])
	if i == 3:
		ax[8+i].set_xlabel('X(km)')
	else:
		ax[8+i].set_xticklabels([])
	ax[i+8].set_yticklabels([])
	


plt.savefig(fig_dir + '2wave_traj_euler.png')

