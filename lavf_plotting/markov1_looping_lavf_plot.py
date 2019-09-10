#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit,fsolve,root
import sympy as sp
from scipy.interpolate import interp1d



import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import square_grid_plot, spd_grid_plot

# INPUT PARAMETERS #

regime = 1

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/LAVF/UNIFORM/' % regime
nbins = 10
ntau = 100

# READ SPD #

R = np.zeros((2,nbins,2,2,ntau))

file_name = home_dir + '/STATS/LAVF/full_new_LAVF.nc'
R_data = Dataset(file_name,'r')
tmp = np.transpose(R_data.variables['LAVF'][:])
R[0,:,:,:,:] = tmp

file_name = home_dir + '/STATS/LAVF/looping_LAVF.nc'
R_data = Dataset(file_name,'r')
tmp = np.transpose(R_data.variables['LAVF'][:])
R[1,:,:,:,:] = tmp

			
# PLOT R

t = np.arange(0,ntau,1)

nrows = 2
ncols = 2
right_space = 0.02
left_space = 0.12
top_space = 0.05
bottom_space = 0.08
hor_space = 0.08
ver_space = 0.08
fig_width = 8.27

for b in range(nbins):
	ax = spd_grid_plot(left_space,right_space,bottom_space,top_space,ver_space,hor_space)
	k=0
	for i in range(ncols):
		for j in range(nrows):
			ax[k].plot(t,R[0,b,i,j,:],'b-',label = 'Full')
			ax[k].plot(t,R[1,b,i,j,:],'g--',label = 'Markov1 Looping')
	
			ax[k].ticklabel_format(style = 'sci',axis = 'y',scilimits=(0,0))
			ax[k].grid()
			if i == 0:
				ax[k].set_ylim([-.2,1])
			else:
				ax[k].set_ylim([-.75,1])
			k+=1
	if (b == 0 or b == 5):
		ax[0].legend(loc = 'lower left')
	
	ax[1].set_xlabel('Time Lag (days)')
	ax[3].set_xlabel('Time Lag (days)')
	ax[0].set_ylabel('Top Layer')
	ax[1].set_ylabel('Bottom Layer') 
	ax[0].set_title('R$_x$')
	ax[2].set_title('R$_y$')
	fig_name = fig_dir + 'markov1_R_bin%i' % (b+1)
	plt.savefig(fig_name)
	plt.close()





	
	




