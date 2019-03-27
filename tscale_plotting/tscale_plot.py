#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import square_grid_plot

# INPUT PARAMETERS #

regime = 2

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i/' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/TSCALE/' % regime
nbins = 10
nt = 1000

# READ TIME-AVERAGED STREAM FUNCTION

file_name = home_dir + 'QG/QG_ave.nc'
psi_data = Dataset(file_name,'r')
psi_ave = np.transpose(psi_data.variables['Time-averaged Stream Function'][:])

# READ ALPHA

tscale = np.zeros((3,nbins,2))

for i in range(3):
	file_name = home_dir + 'STATS/TSCALE/full_MEAN_TSCALE.nc'
	tscale_data = Dataset(file_name,'r')
	tmp = np.transpose(tscale_data.variables['Time Scale'][:])
	tscale[0,:,:] = tmp
	
	file_name = home_dir + 'STATS/TSCALE/eddy_MEAN_TSCALE.nc'
	tscale_data = Dataset(file_name,'r')
	tmp = np.transpose(tscale_data.variables['Time Scale'][:])
	tscale[1,:,:] = tmp
	
	file_name = home_dir + 'STATS/TSCALE/pseudo_MEAN_TSCALE.nc'
	tscale_data = Dataset(file_name,'r')
	tmp = np.transpose(tscale_data.variables['Time Scale'][:])
	tscale[2,:,:] = tmp

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
hor_space = .03
ver_space = 0.
top_space = .15
bottom_space = .15
left_space = .08
right_space = .02
fig_width = 8.27

ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

for i in range(ncols):
	ax[i].pcolor(xx,yy,psi_ave[:,:,i],alpha=0.5,cmap = cm.gray)
	ax_alpha = ax[i].twiny()
	ax_alpha.plot(tscale[0,:,i],bin_centres,'b-',label='Full')
	ax_alpha.plot(tscale[1,:,i],bin_centres,'g--',label='Eddy')
	ax_alpha.plot(tscale[2,:,i],bin_centres,'r-.',label='FFE')
	ax_alpha.grid()
	if (i == 0):
		#ax_alpha.set_title('Top Layer')
		ax[i].set_ylabel('Y (km)')
	else:
		#ax_alpha.set_title('Bottom Layer')
		ax[i].set_yticklabels([])
	ax_alpha.set_xlim(0,1000)
	ax[i].set_ylim(0,520)
	ax_alpha.set_xlabel('Time Scale (Days)')
	ax[i].set_xlabel('X (km)')
	
		
	ax_alpha.legend()

fig_name = fig_dir + 'tscale_uniform'
plt.savefig(fig_name)	

	



