## CALCULATE THE EDDY KINETIC ENERGY ##

#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import square_grid_plot,a4_plot,spd_grid_plot
from calculate_var import velocity_from_psi

plt.rcParams.update({'font.size': 8})

# INPUT PARAMETERS #

regime = 1

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i/' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/QG/ENERGY/' % regime

# READ VELOCITY VARIANCE #

data_dir = home_dir + 'STATS/SIGMA/'
sigma_file = data_dir + 'velocity_variance.nc'
sigma_data = Dataset(sigma_file,'r')
sigma = sigma_data.variables['Velocity Variance'][:]

# READ TRAJECTORIES #

traj_file = home_dir + 'TRAJ/UNIFORM_BINS/pseudo_uniform_bins_trajectories.nc'
traj_data = Dataset(traj_file,'r')
traj = 
