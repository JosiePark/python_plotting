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
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/THETA/' % regime
nbins = 10
nt = 1000

# READ OMEGA DISTRIBUTION #

omega_dir = home_dir + 'STATS/THETA/pseudo_omega_distribution.nc'
omega_data = Dataset(omega_dir,'r')
omega = omega_data.variables['Omega'][:]

hor_space = 0.06
ver_space = 0.03
nrows = 5
ncols = 2
left_space = 0.08
right_space = 0.02
bottom_space = 0.05
top_space = 0.02

ax = a4_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,ver_space,hor_space)

print(min(omega[:,0,2]))
print(max(omega[:,0,2]))

k = 0

for i in range(ncols):
	for j in range(nrows):
		
		ax[k].hist(omega[:,0,k],100,color='k')
		#ax[k].ticklabel_format(style = 'sci',axis = 'x',scilimits=(0,2))
		ax[k].set_ylim([0,3000])
		ax[k].set_xlim([-0.08,0.08])
		if (i == 0):
			ax[k].set_ylabel('Number of Particles')
		if (j == 4) or (j == 9):
			ax[k].set_xlabel('$\Omega$ (days$^{-1}$)')
		ax[k].text(.9,.9,'Bin %i' % (k+1),horizontalalignment = 'center',verticalalignment = 'center', transform = ax[k].transAxes)
		ax[k].grid()
		k+=1


plt.savefig(fig_dir + 'omega_toplayer')
plt.show()
plt.close()

ax = a4_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,ver_space,hor_space)
k = 0

for i in range(ncols):
	for j in range(nrows):
		
		ax[k].hist(omega[:,1,k],100,color='k')
		#ax[k].ticklabel_format(style = 'sci',axis = 'x',scilimits=(0,2))
		ax[k].set_ylim([0,3000])
		ax[k].set_xlim([-0.08,0.08])
		if (i == 0):
			ax[k].set_ylabel('Number of Particles')
		if (j == 4) or (j == 9):
			ax[k].set_xlabel('$\Omega$ (days$^{-1}$)')
		ax[k].text(.9,.9,'Bin %i' % (k+1),horizontalalignment = 'center',verticalalignment = 'center', transform = ax[k].transAxes)
		ax[k].grid()
		k+=1

plt.savefig(fig_dir + 'omega_bottomlayer')
plt.show()
