#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

# INPUT PARAMETERS #

regime = 1

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i/' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/ALPHA/' % regime
nbins = 10
nt = 1000

# READ TIME-AVERAGED STREAM FUNCTION

file_name = home_dir + 'QG/QG_ave.nc'
psi_data = Dataset(file_name,'r')
psi_ave = np.transpose(psi_data.variables['Time-averaged Stream Function'][:])

# READ ALPHA

alpha = np.zeros((3,2,nbins,2))

for i in range(3):
	file_name = home_dir + 'STATS/ALPHA/full_MEAN_ALPHA.nc'
	alpha_data = Dataset(file_name,'r')
	tmp = np.transpose(alpha_data.variables['Alpha'][:])
	alpha[0,:,:,:] = tmp
	
	file_name = home_dir + 'STATS/ALPHA/eddy_MEAN_ALPHA.nc'
	alpha_data = Dataset(file_name,'r')
	tmp = np.transpose(alpha_data.variables['Alpha'][:])
	alpha[1,:,:,:] = tmp
	
	file_name = home_dir + 'STATS/ALPHA/pseudo_MEAN_ALPHA.nc'
	alpha_data = Dataset(file_name,'r')
	tmp = np.transpose(alpha_data.variables['Alpha'][:])
	alpha[2,:,:,:] = tmp

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
		ax_alpha = ax[2*i+j].twiny()
		ax_alpha.plot(alpha[0,j,:,i],bin_centres,'b-',label='Full')
		ax_alpha.plot(alpha[1,j,:,i],bin_centres,'g--',label='Eddy')
		ax_alpha.plot(alpha[2,j,:,i],bin_centres,'r-.',label='FFE')
		ax_alpha.grid()
		if (i == 0):
			if (j == 0):
				ax_alpha.set_xlabel(r'$\alpha_x$')
			else:
				ax_alpha.set_xlabel(r'$\alpha_y$')
		else:
			ax[2*i+j].set_xlabel('X (km)')
		ax_alpha.set_xlim(0,3)
		ax[2*i+j].set_ylim(0,520)
		if (j==0):
			ax[2*i+j].set_ylabel('Y (km)')
		
		ax_alpha.legend()

fig_name = fig_dir + 'alpha_uniform'
plt.savefig(fig_name)

	



