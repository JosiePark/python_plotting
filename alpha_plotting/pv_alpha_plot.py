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

alpha = np.zeros((2,nbins))

for i in range(2):
	file_name = home_dir + 'STATS/ALPHA/full_MEAN_ALPHA.nc'
	alpha_data = Dataset(file_name,'r')
	tmp = np.transpose(alpha_data.variables['Alpha'][0,:,1])
	alpha[0,:] = tmp
	
	file_name = home_dir + 'STATS/ALPHA/test_full_new_PVDISP_MEAN_ALPHA.nc'
	alpha_data = Dataset(file_name,'r')
	tmp = np.transpose(alpha_data.variables['Alpha'][:])
	alpha[1,:] = tmp
	
# READ PV BIN WIDTHS

file_name = home_dir + 'TRAJ/PV_BINS/bin_width.nc'
width_data = Dataset(file_name,'r')
pvbin_width = np.transpose(width_data.variables['Bin Width'][:])
	

# PLOT ALPHA SUPERIMPOSED ON THE TIME-AVERAGED STREAM FUNCTION

xmin,xmax = 0,520

xx,yy = np.mgrid[xmin:xmax:512j,xmin:xmax:512j]

bin_width = 512./10.
bin_centres = []
bin_centres.append(bin_width/2.)
for b in range(nbins-1):
	bin_centres.append(bin_centres[b]+bin_width)
	
# CALCULATE PV BINS

pvbin_centres = []
pvbin_centres.append(pvbin_width[0]/2.)
for b in range(nbins-1):
	pvbin_centres.append(pvbin_centres[b]+(pvbin_width[b] + pvbin_width[b+1])/2.)

print(pvbin_centres)

nrows = 1
ncols = 1	
hor_space = .03
ver_space = 0.04
top_space = .05
bottom_space = .051
left_space = .08
right_space = .02
fig_width = 8.27



k = 0
ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

for i in range(ncols):
	for j in range(nrows):
		
		ax[k].pcolor(xx,yy,psi_ave[:,:,j],alpha=0.5,cmap = cm.gray)
		ax_K = ax[k].twiny()
		ax_K.plot(alpha[0,:],bin_centres,'b-',label='Full')
		ax_K.plot(alpha[1,:],pvbin_centres,'g--',label='PV Mapped')
		ax_K.set_xlabel('alpha$_y$')
		ax[k].set_xlabel('X (km)')
		ax_K.set_xlim(0,3)
		ax[k].set_ylim(0,520)
		ax[k].set_xlim(0,520)
		ax[k].set_ylabel('Y (km)')

		ax_K.legend(loc = 'upper left')
		ax_K.grid()
		k+=1

fig_name = fig_dir + 'alpha_PV'
plt.savefig(fig_name)
print(alpha[1,:])



	



