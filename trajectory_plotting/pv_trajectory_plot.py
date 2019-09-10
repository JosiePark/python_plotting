#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import square_grid_plot,a4_plot,spd_grid_plot

plt.rcParams.update({'font.size': 8})

# NAME TRAJECTORY FILE 

regime = 2
ii = 512
basinscale = 520.
scale = basinscale/float(ii)

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i/' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/TRAJ/' % regime

traj_file = home_dir + 'TRAJ/PV_BINS/full_PV_bins_trajectories.nc'
traj_data = Dataset(traj_file,'r')

#time = traj_data.variables['Time'][:]

# FIND LOOPING TRAJETORY

bin_no = [4,9]
rel_no = 0

print(traj_data)

traj = np.transpose(traj_data.variables['Trajectories'][:,0,rel_no,:,:,bin_no])
coord = np.transpose(traj_data.variables['Domain Coordinate'][:,0,rel_no,:,:,bin_no])

traj = scale*(traj + ii*coord)

print(traj.shape)

nrows = 1
ncols = 2
left_space = 0.08
right_space = 0.02
bottom_space = 0.12
top_space = 0.1
hor_space = 0.08
ver_space = 0.05
fig_width = 8.27

ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

for k in range(2):
	ax[k].plot(np.transpose(traj[k,0,:4,:]),np.transpose(traj[k,1,:4,:]))
	ax[k].set_xlabel('X (km)')
	ax[k].set_ylabel('Y (km)')
	ax[k].set_title('Bin %i' % (bin_no[k] + 1))

plt.savefig(fig_dir + 'pv_traj')
plt.show()
