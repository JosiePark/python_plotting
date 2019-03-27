#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import square_grid_plot,a4_plot

# INPUT PARAMETERS #

regime = 1

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/SPD/DIFFUSION_MODELS/' % regime
nbins = 10
nt = 1000

# READ SPD #

SPD = np.zeros((3,nbins,nt))

for b in range(nbins):
	file_name = home_dir + '/STATS/SPD/FULL/full_uniform_SPD_bin%i.nc' % (b+1)
	spd_data = Dataset(file_name,'r')
	tmp = np.transpose(spd_data.variables['SPD'][:])
	SPD[0,b,:] = np.mean(tmp[1,0,:,:nt],0)
	
	file_name = home_dir + '/STATS/SPD/DIFFUSION/FULL/full_diffusion_SPD_bin%i.nc' % (b+1)
	spd_data = Dataset(file_name,'r')
	tmp = np.transpose(spd_data.variables['SPD'][:])
	SPD[1,b,:] = tmp[1,0,:nt]

	file_name = home_dir + '/STATS/SPD/DIFFUSION/FULL_PV/full_pv_diffusion_SPD_bin%i.nc' % (b+1)
	spd_data = Dataset(file_name,'r')
	tmp = np.transpose(spd_data.variables['SPD'][:])
	SPD[2,b,:] = tmp[:nt]
			
# PLOT SPD

t = np.arange(0,1000,1)

#for b in range(nbins):
#	fig = plt.figure(constrained_layout = True)
#	gs = gridspec.GridSpec(ncols=1,nrows=1,figure=fig)
#	ax = []
#	ax.append(plt.subplot(gs[0]))
#	ax[0].plot(t,SPD[0,b,:],label = 'Full')
#	ax[0].plot(t,SPD[1,b,:],label = 'Diffusion')
#	ax[0].plot(t,SPD[2,b,:],label = 'PV Mapped Diffusion')
#	ax[0].legend()
#	ax[0].ticklabel_format(style = 'sci',axis = 'y',scilimits=(0,0))
#	ax[0].set_xlabel('Time (days)')
#	ax[0].set_ylabel('D$_y$ (km s $^{-1}$)')
#	fig_name = fig_dir + 'diffusion_SPD_bin%i' % (b+1)
#	plt.savefig(fig_name)

hor_space = 0.0
ver_space = 0.02
nrows = 5
ncols = 2
left_space = 0.1
right_space = 0.01
bottom_space = 0.05
top_space = 0.02

ax = a4_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,ver_space,hor_space)

k = 0

for i in range(ncols):
	for j in range(nrows):
		ax[k].plot(t,SPD[0,k,:],label = 'Full')
		ax[k].plot(t,SPD[1,k,:],label = 'Diffusion')
		ax[k].plot(t,SPD[2,k,:],label = 'PV Mapped Diffusion')
		ax[k].legend()
		ax[k].ticklabel_format(style = 'sci',axis = 'y',scilimits=(0,0))	
		ax[k].set_ylim([0,2.e4])
		ax[k].text(.9,.1,'Bin %i' % (k+1),horizontalalignment = 'center',verticalalignment = 'center', transform = ax[k].transAxes)
		if i == 0:
			ax[k].set_ylabel('D$_y$ (km s $^{-1}$)')
		else:
			ax[k].set_yticklabels([])
		if j == nrows-1:
			ax[k].set_xlabel('Time (days)')
		else:
			ax[k].set_xticklabels([])
		k+=1

fig_name = fig_dir + 'diffusion_SPD'
plt.savefig(fig_name)
	
	




