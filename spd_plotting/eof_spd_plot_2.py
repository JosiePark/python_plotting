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

SPD = np.zeros((2,2,nbins,nrel,2,nt))

if regime == 1:

	for b in range(nbins):
		file_name = home_dir + '/STATS/SPD/PSEUDO/pseudo_uniform_SPD_bin%i.nc' % (b+1)
		spd_data = Dataset(file_name,'r')
		print(spd_data)
		tmp = np.transpose(spd_data.variables['SPD'][:nt])
		SPD[0,:,b,:,:,:] = tmp
	
		file_name = home_dir + '/STATS/SPD/EOF_1-10/EOF_1-10_full_uniform_SPD_bin%i.nc' % (b+1)
		spd_data = Dataset(file_name,'r')
		print(spd_data)
		tmp = np.transpose(spd_data.variables['SPD'][:nt])
		SPD[1,:,b,:,:,:] = tmp

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



ave_SPD = np.mean(SPD,3)
			
# PLOT SPD

t = np.arange(0,nt,1)

right_space = 0.02
left_space = 0.08
top_space = 0.05
bottom_space = 0.08
hor_space = 0.06
ver_space = 0.06
nrows = 2
ncols = 2
fig_width = 8.27


if regime == 1:

	for b in range(nbins):
		ax = spd_grid_plot(left_space,right_space,bottom_space,top_space,hor_space,ver_space)
		k=0
		for i in range(ncols):
			for j in range(nrows):
				ax[k].plot(t,ave_SPD[0,i,b,j,:],label = 'FFE')
				ax[k].plot(t,ave_SPD[1,i,b,j,:],label = 'EOFs 1-10')

				if (i == 0 and j == 0):
					ax[k].legend(loc = 'upper left')
				ax[k].ticklabel_format(style = 'sci',axis = 'y',scilimits=(0,0))
				ax[k].grid()
			
				if (j==1):
					ax[k].set_xlabel('Time (days)')
				else:
					if (i==0):
						ax[k].set_title('D$_x$ (km s $^{-1}$)')
					else:
						ax[k].set_title('D$_y$ (km s $^{-1}$)')
				
				if (i==0):
					if (j == 1):
						ax[k].set_ylabel('Bottom Layer')
					else:
						ax[k].set_ylabel('Top Layer')
				k+=1
		fig_name = fig_dir + 'SPD_1-10_bin%i' % (b+1)
		plt.savefig(fig_name)

else:

	for b in range(nbins):
		ax = spd_grid_plot(left_space,right_space,bottom_space,top_space,hor_space,ver_space)
		k=0
		for i in range(ncols):
			for j in range(nrows):
				ax[k].plot(t,ave_SPD[0,i,b,j,:],label = 'FFE')
				ax[k].plot(t,ave_SPD[4,i,b,j,:],label = 'rossbyHalf')
				if (i == 0 and j == 0):
					ax[k].legend(loc = 'upper left')
				ax[k].ticklabel_format(style = 'sci',axis = 'y',scilimits=(0,0))
				ax[k].grid()
			
				if (j==1):
					ax[k].set_xlabel('Time (days)')
				else:
					if (i==0):
						ax[k].set_title('D$_x$ (km s $^{-1}$)')
					else:
						ax[k].set_title('D$_y$ (km s $^{-1}$)')
				
				if (i==0):
					if (j == 1):
						ax[k].set_ylabel('Bottom Layer')
					else:
						ax[k].set_ylabel('Top Layer')
				k+=1
		fig_name = fig_dir + 'SPD_bin%i' % (b+1)
		plt.savefig(fig_name)


	
	




