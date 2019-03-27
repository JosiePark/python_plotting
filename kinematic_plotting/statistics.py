#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import square_grid_plot, spd_grid_plot

plt.rcParams.update({'font.size': 8})

# INPUT PARAMETERS #

regime = 1

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i' % regime
#fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/SPD/KINEMATIC/' % regime
nbins = 10
nt = 1000
nrel = 9

# read files

file_name = home_dir + '/TRAJ/KINEMATIC/traj_top.nc'
traj_data = Dataset(file_name,'r')

file_name = home_dir + '/TRAJ/KINEMATIC/traj_zonalRossby_top.nc'
traj_data = Dataset(file_name,'r')

file_name = home_dir + '/TRAJ/KINEMATIC/traj_top_stationary.nc'
traj_data = Dataset(file_name,'r')

# plot trajectories

# calculate variance, skewness and kurtosis for each bin


