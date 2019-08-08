#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation

import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')

regime = 1
ii = 512
jj = 512

# read data

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/TRAJ/' % regime

traj_file = home_dir + '/TRAJ/DIFFUSION/markov1.nc'
traj_data = Dataset(traj_file,'r')
traj = np.transpose(traj_data.variables['Trajectories'][:])

print(traj_data)
print(traj.shape)

[ndims,npoints,nbins,nlayers,nt] = traj.shape
traj = np.reshape(traj,(ndims,npoints*nbins,nlayers,nt))

# produce movie of scatter plots

def _update_plot(i):
	scat.set_offsets(traj[:,:,0,i])
	print('Frames: %d' %i)
	return scat,

plt.ion()	
plt.figure()
for t in range(nt):
	print('time t = ',nt)
	plt.clf()
	plt.scatter(traj[0,:,0,t] % 512,traj[1,:,0,t] % 512)
	plt.axis((0,512,0,512))
	plt.pause(.01)




