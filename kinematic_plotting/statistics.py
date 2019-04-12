#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import kurtosis,skew


import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import square_grid_plot, spd_grid_plot

plt.rcParams.update({'font.size': 8})

# INPUT PARAMETERS #

regime = 1

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/KINEMATIC/' % regime
nbins = 10
nt = 1000
nrel = 9
npoints = 2500

layer = 2

# read files

traj = np.zeros((4,nbins,npoints,2,nt))

if layer == 1:

	file_name = home_dir + '/TRAJ/KINEMATIC/traj_top.nc'
	traj_data = Dataset(file_name,'r')
	traj[0,:,:,:,:] = np.transpose(traj_data.variables['Trajectories'][:])


	file_name = home_dir + '/TRAJ/KINEMATIC/traj_top_c1.nc'
	traj_data = Dataset(file_name,'r')
	traj[1,:,:,:,:] = np.transpose(traj_data.variables['Trajectories'][:])

	file_name = home_dir + '/TRAJ/KINEMATIC/traj_top_smallc.nc'
	traj_data = Dataset(file_name,'r')
	traj[2,:,:,:,:] = np.transpose(traj_data.variables['Trajectories'][:])
	
else:

	file_name = home_dir + '/TRAJ/KINEMATIC/traj_bottom.nc'
	traj_data = Dataset(file_name,'r')
	traj[0,:,:,:,:] = np.transpose(traj_data.variables['Trajectories'][:])


	file_name = home_dir + '/TRAJ/KINEMATIC/traj_bottom_c1.nc'
	traj_data = Dataset(file_name,'r')
	traj[1,:,:,:,:] = np.transpose(traj_data.variables['Trajectories'][:])

	file_name = home_dir + '/TRAJ/KINEMATIC/traj_bottom_smallc.nc'
	traj_data = Dataset(file_name,'r')
	traj[2,:,:,:,:] = np.transpose(traj_data.variables['Trajectories'][:])


# READ PSEUDO trajectories

file_name = home_dir +'/TRAJ/UNIFORM_BINS/pseudo_uniform_bins_trajectories.nc'
traj_data = Dataset(file_name,'r')
tmp = np.transpose(traj_data.variables['Trajectories'][:,layer-1,:,:300,:,:])
coord = np.transpose(traj_data.variables['Domain Coordinate'][:,layer-1,:,:300,:,:])
tmp = np.reshape(tmp,(nbins,2,300*9,nt))
tmp = np.swapaxes(tmp,1,2)

tmp_coord = np.reshape(tmp,(nbins,2,300*9,nt))
tmp_coord = np.swapaxes(tmp_coord,1,2)

traj[3,:,:,:,:]  = tmp[:,:npoints,:,:] + 512.*tmp_coord[:,:npoints,:,:]
print(traj[3,0,2,:,:])

traj_plot = np.reshape(traj[:,:,:50,:,:],(4,50*nbins,2,nt))


# plot trajectories

nrows = 1
ncols = 3
left_space = 0.1
right_space = 0.01
bottom_space = 0.15
top_space = 0.1
ver_space = 0.
hor_space = 0.
fig_width = 8.27

ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

k = 0

print(traj_plot[1,:,0,:].shape)

for c in range(ncols):
	for r in range(nrows):
		ax[k].plot(np.transpose(traj_plot[c,:,0,:]),np.transpose(traj_plot[c,:,1,:]))
		ax[k].set_xlabel('X(km)')
		k+=1
		

ax[0].set_title('rossbyHalf')
ax[1].set_title('c = 1')
ax[2].set_title('c = 0.01')
#ax[3].set_title('FFE')
ax[0].set_ylabel('Y (km)')	
ax[1].set_yticklabels([])
ax[2].set_yticklabels([])
#ax[3].set_yticklabels([])

if layer == 1:
	fig_name = fig_dir + 'top_trajectories_c'
else:
	fig_name = fig_dir + 'bottom_trajectories_c'
plt.savefig(fig_name)
plt.show()
plt.close()

############################################################
# calculate mean,variance, skewness and kurtosis for each bin

x_dispersion = traj[:,:,:,0,-1] - traj[:,:,:,0,0]
moment = np.zeros((4,4,nbins))

# mean

moment[0,:,:] = np.mean(x_dispersion,axis = 2)

# variance

moment[1,:,:] = np.var(x_dispersion,axis = 2)

# skewness

moment[2,:,:] = skew(x_dispersion,axis = 2)

# kurtosis

moment[3,:,:] = kurtosis(x_dispersion,axis = 2)

nrows = 4
ncols = 3
left_space = 0.05
right_space = 0.01
bottom_space = 0.05
top_space = 0.03
ver_space = 0.05
hor_space = 0.
fig_width = 8.27

ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)
k = 0
for c in range(ncols):
	for r in range(nrows):
		ax[k].plot(np.transpose(moment[r,c,:]),np.linspace(1,nbins,nbins))
		if r == 0:
			ax[k].set_xlabel('Mean')
		elif r == 1:
			ax[k].set_xlabel('Variance')
		elif r == 2:
			ax[k].set_xlabel('Skewness')
		elif r == 3:
			ax[k].set_xlabel('Kurtosis')
		if c != 0:
			ax[k].set_yticklabels([])
		else:
			ax[k].set_ylabel('Bin')
		if r == 0:
			if c == 0:
				ax[k].set_title('rossbyHalf')
			elif c == 1:
				ax[k].set_title('c = 1')
			elif c == 2:
				ax[k].set_title('c = 0.01')
			elif c == 3:
				ax[k].set_title('FFE')
			
		ax[k].grid()
		k+=1
	
if layer == 1:
	fig_name = fig_dir + 'top_statistics_c'
else:
	fig_name = fig_dir + 'bottom_statistics_c'
plt.savefig(fig_name)	
plt.show()

