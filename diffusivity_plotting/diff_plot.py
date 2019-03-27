#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

# INPUT PARAMETERS #

regime = 1

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i/' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/DIFFUSIVITY/' % regime
nbins = 10
nt = 1000

# READ TIME-AVERAGED STREAM FUNCTION

file_name = home_dir + 'QG/QG_ave.nc'
psi_data = Dataset(file_name,'r')
psi_ave = np.transpose(psi_data.variables['Time-averaged Stream Function'][:])

# READ DIFFUSIVITY

K = np.zeros((3,2,nbins,2))

for i in range(3):
	file_name = home_dir + 'STATS/DIFFUSIVITY/full_MEAN_DIFF.nc'
	diff_data = Dataset(file_name,'r')
	tmp = np.transpose(diff_data.variables['Diffusivity'][:])
	K[0,:,:,:] = tmp
	
	file_name = home_dir + 'STATS/DIFFUSIVITY/eddy_MEAN_DIFF.nc'
	diff_data = Dataset(file_name,'r')
	tmp = np.transpose(diff_data.variables['Diffusivity'][:])
	K[1,:,:,:] = tmp
	
	file_name = home_dir + 'STATS/DIFFUSIVITY/pseudo_MEAN_DIFF.nc'
	diff_data = Dataset(file_name,'r')
	tmp = np.transpose(diff_data.variables['Diffusivity'][:])
	K[2,:,:,:] = tmp

# PLOT ALPHA SUPERIMPOSED ON THE TIME-AVERAGED STREAM FUNCTION

xmin,xmax = 0,520

xx,yy = np.mgrid[xmin:xmax:512j,xmin:xmax:512j]

bin_width = 512./10.
bin_centres = []
bin_centres.append(bin_width/2.)
for b in range(nbins-1):
	bin_centres.append(bin_centres[b]+bin_width)


fig = plt.figure(constrained_layout = True)
gs = gridspec.GridSpec(ncols = 2,nrows=2,figure=fig)
ax = []
ax.append(plt.subplot(gs[0]))
ax.append(plt.subplot(gs[1]))
ax.append(plt.subplot(gs[2]))
ax.append(plt.subplot(gs[3]))
for i in range(2):
	for j in range(2):
		ax[2*i+j].pcolor(xx,yy,psi_ave[:,:,j],alpha=0.5,cmap = cm.gray)
		ax_K = ax[2*i+j].twiny()
		ax_K.plot(K[0,j,:,i],bin_centres,label='Full')
		ax_K.plot(K[1,j,:,i],bin_centres,label='Eddy')
		ax_K.plot(K[2,j,:,i],bin_centres,label='FFE')
		if (i == 0):
			if (j == 0):
				ax_K.set_xlabel(r'$\K_x$')
			else:
				ax_K.set_xlabel(r'$\K_y$')
		else:
			ax[2*i+j].set_xlabel('X (km)')
		#ax_K.set_xlim(0,3)
		ax[2*i+j].set_ylim(0,520)
		if (j==0):
			ax[2*i+j].set_ylabel('Y (km)')
		
		ax_alpha.legend()

fig_name = fig_dir + 'K_uniform'
plt.savefig(fig_name)

	



