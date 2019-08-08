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

i_diff = 1
i_theta = 1
i_sigma = 1

# INPUT PARAMETERS #

regime = 1

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i/' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/MARKOV1/' % regime
nbins = 10

# READ PARAMETERS #

theta_dir = home_dir + 'STATS/THETA/'
theta_file = theta_dir + 'pseudo_theta_osc.nc'
theta_data = Dataset(theta_file,'r')
theta = theta_data.variables['Theta'][:]

K_dir = home_dir + 'STATS/DIFFUSIVITY/'
K_file = K_dir + 'pseudo_MEAN_DIFF.nc'
K_data = Dataset(K_file,'r')
K = np.transpose(K_data.variables['Diffusivity'][:])
K = np.swapaxes(K,0,1)

sigma_dir = home_dir + 'STATS/SIGMA/'
sigma_file = sigma_dir + 'velocity_variance.nc'
sigma_data = Dataset(sigma_file,'r')
sigma = sigma_data.variables['Velocity Variance'][:]
#print(sigma[0,:,:,0,0])
plt.contourf(sigma[0,:,:,0,0])
plt.colorbar()
plt.show()

# SCALING 

basinscale = 520.e5
ii = 512
scale = basinscale/float(ii)
uscale = 1.
tscale = scale/uscale

scale_new = scale/(10**5)
tscale_new = tscale/86400.

theta = theta*86400./tscale
K = K*(tscale/86400.)/((scale/(10**5))**2)


# AVERAGE SIGMA ACROSS THE BINS

bin_width = int(512./nbins)

sigma_ave = np.zeros((nbins,2,2)) # bin, dimension, layer
for b in range(nbins):
	y_start = b*bin_width
	y_end = (b+1)*bin_width
	print(y_start,y_end)
	if (y_end >= 512):
		y_end = 511
	for d in range(2):
		for l in range(2):
			sigma_ave[b,d,l] = np.mean(sigma[l,y_start:y_end,:,d,d])
	

# CALCULATE AND PLOT MISSING PARAMETER #

if i_diff == 1:
	K_new = np.multiply(sigma_ave,theta)
	plt.plot(K_new[:,1,0],label='derived using K=sigma*T')
	plt.plot(K[:,1,0],label = 'known from either D or R')
	ratio = np.divide(K_new,K)
	plt.plot(ratio[:,0,0],label = 'ratio')
	plt.legend()
	plt.title('K')
	plt.savefig(fig_dir + 'K.png')
	plt.show()
	
if i_theta == 1:
	theta_new = np.divide(K,sigma_ave)
	plt.plot(theta_new[:,1,0],label='derived')
	plt.plot(theta[:,1,0],label = 'known')
	ratio = np.divide(theta_new,theta)
	plt.plot(ratio[:,1,0],label = 'ratio')
	plt.legend()
	plt.title('T')
	plt.savefig(fig_dir + 'theta.png')
	plt.show()
	
if i_sigma == 1:
	sigma_new = np.divide(K,theta)
	# plot
	plt.plot(sigma_new[:,1,0],label='derived')
	plt.plot(sigma_ave[:,1,0],label = 'known')
	ratio = np.divide(sigma_new,sigma_ave)
	plt.plot(ratio[:,0,0],label = 'ratio')
	plt.legend()
	plt.title('Sigma')
	plt.savefig(fig_dir + 'sigma.png')
	plt.show()
	
# WRITE NEW DIFFUSIVITY TO FILE

file_name = home_dir + '/STATS/DIFFUSIVITY/derived_diffusivity.nc'
K_data = Dataset(file_name,'w',format='NETCDF4_CLASSIC')
K_data.createDimension('Bin',nbins)
K_data.createDimension('Dimension',2)
K_data.createDimension('Layer',2)
K_var = K_data.createVariable('Diffusivity',np.float64,('Bin','Dimension','Layer',))
K_var[:] = K_new

# READ LAGRANGIAN VELOCITY VARIANCE

file_name = home_dir + '/STATS/SIGMA/lagrangian_sigma.nc'
sigma_data = Dataset(file_name,'r')
sigma_L = sigma_data.variables['Lagrangian Velocity Variance'][:]
sigma_L_new = np.zeros((10,2,2))
sigma_L_new[:,0,:] = sigma_L[:,0,0,:]
sigma_L_new[:,1,:] = sigma_L[:,1,1,:]

K_new = np.multiply(sigma_L_new,theta)
plt.plot(K_new[:,1,0],label='derived using K=sigma_L*T')
plt.plot(K[:,1,0],label = 'known from either D or R')
plt.plot(ratio[:,0,0],label = 'ratio')
plt.legend()
plt.title('K')
plt.savefig(fig_dir + 'K_lag.png')
plt.show()
plt.close()


sigma_L_new = sigma_L_new/scale_new**2

plt.plot(sigma_L_new[:,1,0],label='Lagrangian')
plt.plot(sigma_ave[:,1,0],label = 'Eulerian')
plt.legend()
plt.title('Sigma')
plt.savefig(fig_dir + 'sigma_L.png')
plt.show()
	



	

