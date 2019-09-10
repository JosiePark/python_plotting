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

# INPUT PARAMETERS #

regime = 2

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i/' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/THETA/' % regime
nbins = 10
nt = 1000

# READ TIME-AVERAGED STREAM FUNCTION

file_name = home_dir + 'QG/QG_ave.nc'
psi_data = Dataset(file_name,'r')
psi_ave = np.transpose(psi_data.variables['Time-averaged Stream Function'][:])

# READ ALPHA

file_name = home_dir + 'STATS/THETA/omega_estimate.nc'
alpha_data = Dataset(file_name,'r')
tmp = np.transpose(alpha_data.variables['Omega'][:])
print(tmp.shape)
print(alpha_data)
omega = tmp

basinscale = 520.e5
scale = basinscale/512

	

# PLOT ALPHA SUPERIMPOSED ON THE TIME-AVERAGED STREAM FUNCTION

xmin,xmax = 0,520

xx,yy = np.mgrid[xmin:xmax:512j,xmin:xmax:512j]

bin_width = 512./10.
bin_centres = []
bin_centres.append(bin_width/2.)
for b in range(nbins-1):
	bin_centres.append(bin_centres[b]+bin_width)

nrows = 1
ncols = 2	
hor_space = .04
ver_space = 0.04
top_space = .11
bottom_space = .1
left_space = .08
right_space = .02
fig_width = 8.27



k = 0
ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

for i in range(ncols):
	for j in range(nrows):
		
		ax[k].pcolor(xx,yy,psi_ave[:,:,i],alpha=0.5,cmap = cm.gray)
		ax_K = ax[k].twiny()
		ax_K.plot(omega[:,i],bin_centres,'b-')
		if (j == 0):
			if (i == 0):
				ax_K.set_xlabel('$\Omega$ (days$^{-1}$)')
			else:
				ax_K.set_xlabel('$\Omega$ (days$^{-1}$)')
			ax[k].set_xticklabels([])
		ax[k].set_xlabel('X (km)')
		ax_K.set_xlim(-.035,.02)
		ax[k].set_ylim(0,520)
		ax[k].set_xlim(0,520)
		if (i==0):
			ax[k].set_ylabel('Y (km)')
		else:
			ax[k].set_yticklabels([])
		
		#if (i == 0 and j == 0):
		#	ax_K.legend(loc = 'upper left')
		ax_K.grid()
		k+=1

fig_name = fig_dir + 'omega_uniform'
plt.savefig(fig_name)



	


                             
