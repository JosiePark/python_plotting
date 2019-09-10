#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit,fsolve,root
from scipy.interpolate import interp1d


import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import a4_plot

# INPUT PARAMETERS #
regime = 1

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/LAVF/PV/' % regime
nbins = 10
ntau = 100

# READ SPD #

#file_name = home_dir + '/STATS/LAVF/full_PVDISP_LAVF.nc'
#R_data = Dataset(file_name,'r')
#R1 = np.transpose(R_data.variables['LAVF'][:])

file_name = home_dir + '/STATS/LAVF/full_PVDISP_LAVF.nc'
R_data = Dataset(file_name,'r')
R2 = np.transpose(R_data.variables['LAVF'][:])
#print(R2.shape)

#file_name = home_dir + '/STATS/LAVF/PV_LAVF_2.nc'
#R_data = Dataset(file_name,'r')
#R2 = R_data.variables['LAVF'][:]
			
# PLOT R

t = np.arange(0,ntau,1)

nrows = 1
ncols = 1
right_space = 0.02
left_space = 0.12
top_space = 0.05
bottom_space = 0.08
hor_space = 0.08
ver_space = 0.08
fig_width = 8.27

## FIT R TO AN EXPONENTIAL

def func_x(x,b):
	return np.exp(-b*x)

def func_y(x,a,b):
	return np.exp(-a*x)*(np.cos(b*x))
	
theta = np.zeros((nbins))
omega = np.zeros((nbins))

nrows = 5
ncols = 2
left_space = 0.1
right_space = 0.1
bottom_space = 0.05
top_space = 0.01
right_space = 0.01
hor_space = 0.02
ver_space = 0.0
fig_width = 8.27

ax = a4_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space)

k = 0

t_max = 15
t_int = np.arange(0,t_max,1)

	
for i in range(ncols):
	for j in range(nrows):
		popt, pcov = curve_fit(func_y,t,R2[k,:],maxfev=1500)
		omega[k] = popt[1]
		#popt1,pcov2 = curve_fit(func_y,t_int,R2[k,:t_max],maxfev=1500)
		#f = interp1d(func_y(t_int,*popt1),t_int)
		#theta[k] = f(np.exp(-1))
		theta[k] = 1./popt[0]
		#print(np.exp(-1))
		print(theta[k])
		print('coeffs=',popt)
		ax[k].plot(t,func_y(t,*popt),'r-.',label = "Fitted Curve")
		#omega[k] = popt[1]
		
		if j == nrows-1:
			ax[k].set_xlabel('Time Lag (days)')
		else:
			ax[k].set_xticklabels([])
		if i == 0:
			ax[k].set_ylabel('R$_y$')
		else:
			ax[k].set_yticklabels([])

		#ax[k].plot(t,R1[k,:],'b-',label = "PV Mapped")
		ax[k].plot(t,R2[k,:],'b-',label = "PV Mapped")			
		ax[0].legend(loc = 'upper left')
		ax[k].grid()
		ax[k].text(.5,.05,r'$ T_L =$ %.2f' % theta[k], horizontalalignment='center',verticalalignment='center', transform=ax[k].transAxes)
		ax[k].text(.9,.05,'Bin %i' % (k+1),horizontalalignment = 'center',verticalalignment = 'center', transform = ax[k].transAxes)
		
		ax[k].set_ylim([-1,1])
		
		k+=1

fig_name = fig_dir + 'test_PV_LAVF_osc'
plt.savefig(fig_name)
plt.close()	

	
	
# WRITE THETA TO FILE

file_name = home_dir + '/STATS/THETA/markov1_looping_theta_PV.nc'
theta_data = Dataset(file_name,'w',format='NETCDF4_CLASSIC')
theta_data.createDimension('Bin',nbins)
theta_var = theta_data.createVariable('Theta',np.float64,('Bin'))
theta_var[:] = theta 

file_name = home_dir + '/STATS/THETA/markov1_looping_omega_PV.nc'
omega_data = Dataset(file_name,'w',format = 'NETCDF4_CLASSIC')
omega_data.createDimension('Bin',nbins)
omega_var = omega_data.createVariable('Omega',np.float64,('Bin'))
omega_var[:] = omega





	
	




