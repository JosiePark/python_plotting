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

regime = 1
home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i/' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/TRAJ/' % regime

traj_file = home_dir + 'TRAJ/UNIFORM_BINS/full_uniform_bins_trajectories.nc'
traj_data = Dataset(traj_file,'r')

time = traj_data.variables['Time'][:]

# DETERMINE WHICH TIMES YOU WISH TO READ FROM 

time_read = [0,50,100,150]
[nt,nlayers,nrel,npoints,ndims,nbins] = traj_data.variables['Trajectories'].shape
traj = np.zeros((nbins,ndims,npoints,nlayers,4))

k=0
rel_no = 0
bin_no = 8
for t in time_read:
	tmp = np.transpose(traj_data.variables['Trajectories'][t,:,rel_no,:,:,:])
	traj[:,:,:,:,k] = tmp
	k+=1
	
# SCATTER PLOT

nrows = 2
ncols = 2
left_space = 0.08
right_space = 0.02
bottom_space = 0.05
top_space = 0.05
hor_space = 0.
ver_space = 0.05
fig_width = 8.27

ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)
k=0

read_index = [0,2,1,3]

for k in range(4):
	ax[k].scatter(traj[bin_no,0,:,0,read_index[k]],traj[bin_no,1,:,0,read_index[k]],c = 'k',s = 2)
	ax[k].set_title('Time = %i days' % int(time_read[read_index[k]]))
	ax[k].set_ylim([0,512])
	ax[k].set_xlim([0,512])


ax[0].set_ylabel('Y (km)')
ax[1].set_ylabel('Y (km)')
ax[1].set_xlabel('X (km)')
ax[3].set_xlabel('X (km)')
ax[0].set_xticks([])
ax[2].set_xticks([])
ax[2].set_yticks([])
ax[3].set_yticks([])
fig_name = fig_dir + 'bin%i' % bin_no
plt.savefig(fig_name)

