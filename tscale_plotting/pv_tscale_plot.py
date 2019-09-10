#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import square_grid_plot

plt.rcParams.update({'font.size': 8})

# INPUT PARAMETERS #

regime = 1

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i/' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/TSCALE/' % regime
nbins = 10
nt = 1000

# READ TIME-AVERAGED STREAM FUNCTION

file_name = home_dir + 'QG/QG_ave.nc'
psi_data = Dataset(file_name,'r')
psi_ave = np.transpose(psi_data.variables['Time-averaged Stream Function'][:])

# READ ALPHA

tscale = np.zeros((2,nbins))

for i in range(3):
	file_name = home_dir + 'STATS/TSCALE/full_MEAN_TSCALE.nc'
	tscale_data = Dataset(file_name,'r')
	tmp = np.transpose(tscale_data.variables['Time Scale'][:])
	tscale[0,:] = tmp[:,0]
	
	file_name = home_dir + 'STATS/TSCALE/test_full_PVDISP_MEAN_TSCALE.nc'
	tscale_data = Dataset(file_name,'r')
	tmp = np.transpose(tscale_data.variables['Time Scale'][:])
	tscale[1,:] = tmp

# PLOT ALPHA SUPERIMPOSED ON THE TIME-AVERAGED STREAM FUNCTION

xmin,xmax = 0,520

xx,yy = np.mgrid[xmin:xmax:512j,xmin:xmax:512j]

bin_width = 512./10.
bin_centres = []
bin_centres.append(bin_width/2.)
for b in range(nbins-1):
	bin_centres.append(bin_centres[b]+bin_width)
	
# READ PV BINS #

file_name = home_dir + '/TRAJ/PV_BINS/test_bin_width.nc'
bin_data = Dataset(file_name,'r')
bin_width = np.mean(bin_data.variables['Release Bin Width'][:],0)

bin_boundary = np.zeros((nbins+1))
for b in range(nbins):
	bin_boundary[b+1] = bin_boundary[b] + bin_width[b]
pv_bin_centres = .5*(bin_boundary[:-1]+bin_boundary[1:])
print(pv_bin_centres)

nrows = 1
ncols = 1
hor_space = .03
ver_space = 0.
top_space = .12
bottom_space = .15
left_space = .15
right_space = .04
fig_width = 8.27/2.

ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

for i in range(ncols):
	ax[i].pcolor(xx,yy,psi_ave[:,:,i],alpha=0.5,cmap = cm.gray)
	ax_alpha = ax[i].twiny()
	ax_alpha.plot(tscale[0,:],bin_centres,'b',label='Full')
	ax_alpha.plot(tscale[1,:],pv_bin_centres,'r-.',label='PV mapped')
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
	
	if (i == 0):
		ax_alpha.legend()

fig_name = fig_dir + 'pv_tscale_uniform'
plt.savefig(fig_name)	

	



