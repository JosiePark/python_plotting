#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import square_grid_plot,a4_plot,spd_grid_plot

plt.rcParams.update({'font.size': 8})

# INPUT PARAMETERS #

regime = 1

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/SPD/EOF/' % regime
nbins = 10
nrel = 9
nt = 999

# READ SPD #

SPD = np.zeros((6,2,nbins,nrel,2,nt))

if regime == 1:

	for b in range(nbins):
		file_name = home_dir + '/STATS/SPD/PSEUDO/pseudo_uniform_SPD_bin%i.nc' % (b+1)
		spd_data = Dataset(file_name,'r')
		tmp = np.transpose(spd_data.variables['SPD'][:nt])
		SPD[0,:,b,:,:,:] = tmp
	
		file_name = home_dir + '/STATS/SPD/EOF_MINUS_1-2/EOF_minus_1-2_pseudo_uniform_SPD_bin%i.nc' % (b+1)
		spd_data = Dataset(file_name,'r')
		tmp = np.transpose(spd_data.variables['SPD'][:nt])
		SPD[1,:,b,:,:,:] = tmp

		file_name = home_dir + '/STATS/SPD/EOF_MINUS_3-4/EOF_minus_3-4_pseudo_uniform_SPD_bin%i.nc' % (b+1)
		spd_data = Dataset(file_name,'r')
		tmp = np.transpose(spd_data.variables['SPD'][:nt])
		SPD[2,:,b,:,:,:] = tmp

		file_name = home_dir + '/STATS/SPD/EOF_MINUS_5-6/EOF_minus_5-6_pseudo_uniform_SPD_bin%i.nc' % (b+1)
		spd_data = Dataset(file_name,'r')
		tmp = np.transpose(spd_data.variables['SPD'][:nt])
		SPD[3,:,b,:,:,:] = tmp

		file_name = home_dir + '/STATS/SPD/EOF_MINUS_7-8/EOF_minus_7-8_pseudo_uniform_SPD_bin%i.nc' % (b+1)
		spd_data = Dataset(file_name,'r')
		tmp = np.transpose(spd_data.variables['SPD'][:nt])
		SPD[4,:,b,:,:,:] = tmp

		file_name = home_dir + '/STATS/SPD/EOF_MINUS_9-10/EOF_minus_9-10_pseudo_uniform_SPD_bin%i.nc' % (b+1)
		spd_data = Dataset(file_name,'r')
		tmp = np.transpose(spd_data.variables['SPD'][:nt])
		SPD[5,:,b,:,:,:] = tmp

else:

	for b in range(nbins):
		file_name = home_dir + '/STATS/SPD/PSEUDO/pseudo_uniform_SPD_bin%i.nc' % (b+1)
		spd_data = Dataset(file_name,'r')
		tmp = np.transpose(spd_data.variables['SPD'][:nt])
		SPD[0,:,b,:,:,:] = tmp
	
		file_name = home_dir + '/STATS/SPD/EOF_MINUS_14/EOF_minus_14_pseudo_uniform_SPD_bin%i.nc' % (b+1)
		spd_data = Dataset(file_name,'r')
		tmp = np.transpose(spd_data.variables['SPD'][:nt])
		SPD[1,:,b,:,:,:] = tmp

		file_name = home_dir + '/STATS/SPD/EOF_MINUS_23/EOF_minus_23_pseudo_uniform_SPD_bin%i.nc' % (b+1)
		spd_data = Dataset(file_name,'r')
		tmp = np.transpose(spd_data.variables['SPD'][:nt])
		SPD[2,:,b,:,:,:] = tmp

		file_name = home_dir + '/STATS/SPD/EOF_MINUS_56/EOF_minus_56_pseudo_uniform_SPD_bin%i.nc' % (b+1)
		spd_data = Dataset(file_name,'r')
		tmp = np.transpose(spd_data.variables['SPD'][:nt])
		SPD[3,:,b,:,:,:] = tmp

		file_name = home_dir + '/STATS/SPD/EOF_MINUS_78/EOF_minus_78_pseudo_uniform_SPD_bin%i.nc' % (b+1)
		spd_data = Dataset(file_name,'r')
		tmp = np.transpose(spd_data.variables['SPD'][:nt])
		SPD[4,:,b,:,:,:] = tmp

		file_name = home_dir + '/STATS/SPD/EOF_MINUS_910/EOF_minus_910_pseudo_uniform_SPD_bin%i.nc' % (b+1)
		spd_data = Dataset(file_name,'r')
		tmp = np.transpose(spd_data.variables['SPD'][:nt])
		SPD[5,:,b,:,:,:] = tmp



ave_SPD = np.mean(SPD,3)

# compare contributions across the bins

# READ TIME SCALE #

file_name = home_dir + '/STATS/TSCALE/pseudo_MEAN_TSCALE.nc'
tscale_data = Dataset(file_name,'r')
tscale = np.transpose(tscale_data.variables['Time Scale'][:])

right_space = 0.02
left_space = 0.08
top_space = 0.05
bottom_space = 0.08
hor_space = 0.06
ver_space = 0.06
nrows = 2
ncols = 2
fig_width = 8.27

ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

contribution = np.zeros((5,2,nbins,2))

for l in range(2):
	for b in range(nbins):
		if int(tscale[b,l]) == 1000:
			T = nt-1
		else:
			T = int(tscale[b,l])-1
		for i in range(5):
			contribution[i,:,b,l] = (ave_SPD[i+1,:,b,l,T] - ave_SPD[0,:,b,l,T])/ave_SPD[0,:,b,l,T]



k = 0

for c in range(ncols):
	for r in range(nrows):
		line = []
		if regime == 1:
			ax[k].plot(contribution[0,c,:,r],np.linspace(1,nbins,nbins),label = 'rossbyHalf')
			ax[k].plot(contribution[1,c,:,r],np.linspace(1,nbins,nbins),label = 'rossby')
			ax[k].plot(contribution[2,c,:,r],np.linspace(1,nbins,nbins),label = 'zonal')
			ax[k].plot(contribution[3,c,:,r],np.linspace(1,nbins,nbins),label = 'oscillations')
			ax[k].plot(contribution[4,c,:,r],np.linspace(1,nbins,nbins),label = 'altZonal')
		else:
			ax[k].plot(contribution[0,c,:,r],np.linspace(1,nbins,nbins),label = 'rossbyHalf')
			ax[k].plot(contribution[3,c,:,r],np.linspace(1,nbins,nbins),label = 'rossby')
			ax[k].plot(contribution[1,c,:,r],np.linspace(1,nbins,nbins),label = 'zonal')
			ax[k].plot(contribution[2,c,:,r],np.linspace(1,nbins,nbins),label = 'oscillations')
			ax[k].plot(contribution[4,c,:,r],np.linspace(1,nbins,nbins),label = 'altZonal')

		
		ax[k].set_xlim([-2,2])
		if (r == 1):
			ax[k].set_xlabel('Contribution to SPD')
			if (c == 0):
				ax[k].set_ylabel('Bottom Layer')
		if (r == 0) and (c == 0):
			ax[k].set_ylabel('Top Layer')
			ax[k].set_title('Zonal')
			ax[k].legend(loc = 'upper left')
		if (r == 0) and (c == 1):
			ax[k].set_title('Meridional')
		if (c==1):
			ax[k].set_yticklabels([])
		if (r == 0):
			ax[k].set_xticklabels([])
		ax[k].grid()
		
	

		k+=1

fig_name = fig_dir + 'EOF_contribution'
plt.savefig(fig_name)

plt.show()


	
	




