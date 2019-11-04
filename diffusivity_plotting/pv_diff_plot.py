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
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/DIFFUSIVITY/' % regime
nbins = 10
nt = 1000

# READ TIME-AVERAGED STREAM FUNCTION

file_name = home_dir + 'QG/QG_ave.nc'
psi_data = Dataset(file_name,'r')
psi_ave = np.transpose(psi_data.variables['Time-averaged Stream Function'][:])

# READ DIFFUSIVITY

K = np.zeros((2,nbins))

file_name = home_dir + 'STATS/DIFFUSIVITY/pseudo_new_DIFF.nc'
diff_data = Dataset(file_name,'r')
tmp = np.transpose(diff_data.variables['Diffusivity'][0,:,1])
K[0,:] = tmp

file_name = home_dir + 'STATS/DIFFUSIVITY/test_full_PVDISP_MEAN_DIFF.nc'
diff_data = Dataset(file_name,'r')
tmp = np.transpose(diff_data.variables['Diffusivity'][:])
K[1,:] = tmp

K = K/86400.

K = K*(1000**2)
	

# PLOT ALPHA SUPERIMPOSED ON THE TIME-AVERAGED STREAM FUNCTION

xmin,xmax = 0,520

xx,yy = np.mgrid[xmin:xmax:512j,xmin:xmax:512j]

bin_width = 512./10.
bin_centres = []
bin_centres.append(bin_width/2.)
for b in range(nbins-1):
	bin_centres.append(bin_centres[b]+bin_width)
	
# READ PV BIN WIDTHS # 

file_name = home_dir + '/TRAJ/PV_BINS/test_bin_width.nc'
bin_data = Dataset(file_name,'r')
bin_width = np.mean(bin_data.variables['Release Bin Width'][:],0)

bin_boundary = np.zeros((nbins+1))
for b in range(nbins):
	bin_boundary[b+1] = bin_boundary[b] + bin_width[b]
pv_bin_centres = .5*(bin_boundary[:-1]+bin_boundary[1:])
print(pv_bin_centres)

# CALCULATE BIN CENTRES FOR PV BINS #
	
nrows = 1
ncols = 1
hor_space = .03
ver_space = 0.
top_space = .12
bottom_space = .15
left_space = .15
right_space = .04
fig_width = 8.27/2.

k = 0
ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

for i in range(ncols):
	for j in range(nrows):	
		ax[k].pcolor(xx,yy,psi_ave[:,:,j],alpha=0.5,cmap = cm.gray)
		ax_K = ax[k].twiny()
		ax_K.plot(K[0,:],bin_centres,'b',label='Full')
		ax_K.plot(K[1,:],pv_bin_centres,'r-.',label='PV Mapped')
		if (j == 0):
			if (i == 0):
				ax_K.set_xlabel('K$_y$ (m$^{2}$ s$^{-1}$)')
			else:
				ax_K.set_xlabel('K$_y$ (m$^{2}$ s$^{-1}$)')
			#ax[k].set_xticklabels([])
		ax[k].set_xlabel('X (km)')
		#ax_K.set_xlim(0,3)
		ax[k].set_ylim(0,520)
		ax[k].set_xlim(0,520)
		#ax_K.ticklabel_format(axis='x', style='sci', scilimits=(-2, 2))
		#ax_K.ticklabel_format(style = 'sci',axis = 'x',scilimits=(0,0))
		if (i==0):
			ax[k].set_ylabel('Y (km)')
		else:#
			ax[k].set_yticklabels([])
		
		if (i == 0 and j == 0):
			ax_K.legend(loc = 'upper left')
		ax_K.grid()
		k+=1

fig_name = fig_dir + 'new_K_PV'
plt.savefig(fig_name)

	



