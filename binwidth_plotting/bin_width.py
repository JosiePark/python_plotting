#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import square_grid_plot,a4_plot

# INPUT PARAMETERS #

regime = 2

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/' % regime

# READ BIN WIDTH FILE #

file_name = home_dir + '/TRAJ/PV_BINS/test_bin_width.nc'
bin_data = Dataset(file_name,'r')
bin_width = bin_data.variables['Release Bin Width'][:]
[nrel,nbins] = bin_width.shape

# PLOT BIN WIDTH #

bin_boundary = np.zeros((nrel,nbins+1))

for k in range(nrel):
	for b in range(nbins):
		bin_boundary[k,b+1] = bin_boundary[k,b] + bin_width[k,b]
		
plt.plot(bin_boundary,'k-')
plt.xlabel('Release Number')
plt.ylabel('Y')
plt.title('Bin Width')
fig_name = fig_dir + 'bin_width'
plt.savefig(fig_name)
plt.show()

	




