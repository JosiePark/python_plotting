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

regime = 2

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/LAVF/UNIFORM/' % regime
nbins = 10
ntau = 100

# READ OMEGA # 

file_name = home_dir + '/STATS/THETA/omega_estimate.nc'
omega_data = Dataset(file_name,'r')
omega = omega_data.variables['Omega'][:]
print(omega.shape)

file_name = home_dir + '/STATS/LAVF/pseudo_new_LAVF.nc'
R_data = Dataset(file_name,'r')
R = np.transpose(R_data.variables['LAVF'][:])

			
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

## FIT R TO AN EXPONENTIAL

popt = np.zeros((nbins,2,2))

for b in range(nbins):
	for l in range(2):
		def func_y(x,a):
			return np.exp(-a*x)*(np.cos(omega[l,b]*x))
			
		for d in range(2):
			popt[b,d,l], pcov = curve_fit(func_y,t,R[b,d,l,:])
			
	
theta = 1./popt
	
# WRITE THETA TO FILE

file_name = home_dir + '/STATS/THETA/theta_fitted_looping.nc'
theta_data = Dataset(file_name,'w',format='NETCDF4_CLASSIC')
theta_data.createDimension('Bin',nbins)
theta_data.createDimension('Dimension',2)
theta_data.createDimension('Layer',2)
theta_var = theta_data.createVariable('Theta',np.float64,('Bin','Dimension','Layer',))
theta_var[:] = theta 





	
	




