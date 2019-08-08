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

traj_file = data_dir + '/1/TRAJ/KINEMATIC/traj_top.nc'
traj_data = Dataset(traj_file,'r')
tmp = np.transpose(traj_data.variables['Trajectories'][:])
[nbins,npoints,ndims,nt] = tmp.shape
traj = np.zeros((2,2,nbins,npoints,ndims,nt))
traj[0,0,:,:,:,:] = tmp

traj_file = data_dir + '/1/TRAJ/KINEMATIC/traj_bottom.nc'
traj_data = Dataset(traj_file,'r')
tmp = np.transpose(traj_data.variables['Trajectories'][:])
traj[0,1,:,:,:,:] = tmp

traj_file = data_dir + '/2/TRAJ/KINEMATIC/traj_top.nc'
traj_data = Dataset(traj_file,'r')
tmp = np.transpose(traj_data.variables['Trajectories'][:])
traj[1,0,:,:,:,:] = tmp

traj_file = data_dir + '/2/TRAJ/KINEMATIC/traj_bottom.nc'
traj_data = Dataset(traj_file,'r')
tmp = np.transpose(traj_data.variables['Trajectories'][:])
traj[1,1,:,:,:,:] = tmp

x_dispersion = traj[:,:,:,:,0,-1] - traj[:,:,:,:,0,0]
spread = np.mean(x_dispersion,axis = 3)
print(spread.shape)

# calculate spread in narrower bins

nbins_new = 100
bins = np.linspace(0,512,nbins_new+1)
bins_centre = (bins[:-1]+bins[1:])/2.
spread_new = np.zeros((2,2,nbins_new))
traj_new = np.reshape(traj,[2,2,nbins*npoints,ndims,nt])

bin_index = np.digitize(traj_new[:,:,:,1,0],bins,right = False)
print(bin_index.shape)

for b in range(nbins_new):
	for i in range(2):
		for j in range(2):
			x_dispersion_new = [traj_new[i,j,p,0,-1] - traj_new[i,j,p,0,0] for p in range(len(bin_index[i,j,:])) if bin_index[i,j,p] == b+1]
			spread_new[i,j,b] = np.mean(x_dispersion_new)


traj = np.reshape(traj,[2,2,nbins*npoints,ndims,nt])

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
		ax[k].scatter(np.transpose(traj[i,j,:,0,-1]),np.transpose(traj[i,j,:,1,-1]))
		ax[k].set_xlim([-4000,6000])
		ax[k].grid()
		k+=1	
ax[0].set_title('Coherent Jet')
ax[2].set_title('Latent Jet')
ax[0].set_ylabel('Top Layer')
ax[1].set_ylabel('Bottom Layer')
ax[1].set_xlabel('Stokes Drift')
ax[3].set_xlabel('Stokes Drift')
			
fig_name = fig_dir + 'rossby_half_traj.png'
plt.savefig(fig_name)
plt.show()
plt.close()

ax = square_grid_plot(2,2,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)
k = 0

for i in range(2):
	for j in range(2):
		ax[k].plot(np.transpose(spread_new[i,j,:]),bins_centre,'b-')
		ax[k].set_xlim([-1500,1500])
		ax[k].grid()
		k+=1	
ax[0].set_title('Coherent Jet')
ax[2].set_title('Latent Jet')
ax[0].set_ylabel('Top Layer')
ax[1].set_ylabel('Bottom Layer')
ax[1].set_xlabel('Mean Dispersion')
ax[3].set_xlabel('Mean Dispersion')
			
fig_name = fig_dir + 'rossby_half_dispersion.png'
plt.savefig(fig_name)
plt.show()
plt.close()



