#!/usr/bin/env python

## plot kinematic model trajectories ##

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import square_grid_plot

ii = 512
jj = 512

fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/COMBO/KINEMATIC/'
data_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/'

traj_file = data_dir + '/1/TRAJ/UNIFORM_BINS/eddy_uniform_bins_trajectories.nc'
traj_data = Dataset(traj_file,'r')
tmp = np.transpose(traj_data.variables['Trajectories'][:,:,:,:200,:,:])
[nbins,ndims,npoints,nrel,nlayers,nt] = tmp.shape
tmp = np.swapaxes(tmp,0,4)
tmp = np.swapaxes(tmp,1,4)
print(tmp.shape)
coord = np.transpose(traj_data.variables['Domain Coordinate'][:,:,:,:200,:,:])
coord = np.swapaxes(coord,0,4)
coord = np.swapaxes(coord,1,4)
traj = np.zeros((2,2,nbins,npoints,nrel,ndims,nt))
print(nbins,npoints,nrel,ndims,nt)
traj[0,:,:,:,:,:,:] = tmp + ii*coord

traj_file = data_dir + '/2/TRAJ/UNIFORM_BINS/eddy_uniform_bins_trajectories.nc'
traj_data = Dataset(traj_file,'r')
tmp = np.transpose(traj_data.variables['Trajectories'][:,:,:,:200,:,:])
tmp = np.swapaxes(tmp,0,4)
tmp = np.swapaxes(tmp,1,4)
coord = np.transpose(traj_data.variables['Domain Coordinate'][:,:,:,:200,:,:])
coord = np.swapaxes(coord,0,4)
coord = np.swapaxes(coord,1,4)
traj[1,:,:,:,:,:,:] = tmp + ii*coord

# calculate spread in narrower bins

nbins_new = 100
bins = np.linspace(0,512,nbins_new+1)
bins_centre = (bins[:-1]+bins[1:])/2.
spread_new = np.zeros((2,2,nbins_new))
traj_new = np.reshape(traj,[2,2,nbins*npoints*nrel,ndims,nt])
print(traj_new.shape)

bin_index = np.digitize(traj_new[:,:,:,1,0],bins,right = False)
print(bin_index.shape)

for b in range(nbins_new):
	for i in range(2):
		for j in range(2):
			x_dispersion_new = [traj_new[i,j,p,0,-1] - traj_new[i,j,p,0,0] for p in range(len(bin_index[i,j,:])) if bin_index[i,j,p] == b+1]
			spread_new[i,j,b] = np.mean(x_dispersion_new)
			
left_space = 0.08
right_space = 0.05
bottom_space = 0.08
top_space = 0.05
hor_space = 0.05
ver_space = 0.05
fig_width = 8.27
			
ax = square_grid_plot(2,2,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)
k = 0

for i in range(2):
	for j in range(2):
		ax[k].plot(np.transpose(spread_new[i,j,:]),bins_centre,'b-')
		#ax[k].set_xlim([-1500,1500])
		ax[k].grid()
		k+=1	
ax[0].set_title('Coherent Jet')
ax[2].set_title('Latent Jet')
ax[0].set_ylabel('Top Layer')
ax[1].set_ylabel('Bottom Layer')
ax[1].set_xlabel('Mean Dispersion')
ax[3].set_xlabel('Mean Dispersion')
			
fig_name = fig_dir + 'eddy_dispersion.png'
plt.savefig(fig_name)
plt.show()
plt.close()
