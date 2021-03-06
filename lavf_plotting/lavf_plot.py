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

R = np.zeros((3,nbins,2,2,ntau))

file_name = home_dir + '/STATS/LAVF/full_new_LAVF.nc'
R_data = Dataset(file_name,'r')
tmp = np.transpose(R_data.variables['LAVF'][:])
R[0,:,:,:,:] = tmp

file_name = home_dir + '/STATS/LAVF/eddy_new_LAVF.nc'
R_data = Dataset(file_name,'r')
tmp = np.transpose(R_data.variables['LAVF'][:])
R[1,:,:,:,:] = tmp

file_name = home_dir + '/STATS/LAVF/pseudo_new_LAVF.nc'
R_data = Dataset(file_name,'r')
tmp = np.transpose(R_data.variables['LAVF'][:])
R[2,:,:,:,:] = tmp
			
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
			ax[k].plot(t,R[1,b,i,j,:],'g--',label = 'Eddy')
			ax[k].plot(t,R[2,b,i,j,:],'r-.',label = 'FFE')
	
			ax[k].ticklabel_format(style = 'sci',axis = 'y',scilimits=(0,0))
			ax[k].grid()
			if i == 0:
				ax[k].set_ylim([-.2,1])
			else:
				ax[k].set_ylim([-.75,1])
			k+=1
	if (b == 0 or b == 6):
		ax[0].legend(loc = 'lower left')
	
	ax[1].set_xlabel('Time Lag (days)')
	ax[3].set_xlabel('Time Lag (days)')
	ax[0].set_ylabel('Top Layer')
	ax[1].set_ylabel('Bottom Layer') 
	ax[0].set_title('R$_x$')
	ax[2].set_title('R$_y$')
	fig_name = fig_dir + 'R_bin%i' % (b+1)
	plt.savefig(fig_name)
	plt.close()

## FIT R TO AN EXPONENTIAL

def func_x(x,b):
	return np.exp(-b*x)

def func_y(x,a,b):
	return np.exp(-a*x)*(np.cos(b*x))
	
theta = np.zeros((nbins,2,2))
omega = np.zeros((nbins,2,2))


for b in range(nbins):
	print('bin=',b)
	ax = spd_grid_plot(left_space,right_space,bottom_space,top_space,ver_space,hor_space)
	k=0
	for i in range(ncols):
		for j in range(nrows):
			if i == 0:
				popt, pcov = curve_fit(func_x,t,R[2,b,i,j,:])
				ax[k].plot(t,func_x(t,*popt),'r-.',label = "Fitted Curve")
				theta[b,i,j] = 1./popt[0]
				omega[b,i,j] = popt[0]
			else:
				popt, pcov = curve_fit(func_y,t,R[2,b,i,j,:])
				ax[k].plot(t,func_y(t,*popt),'r-.',label = "Fitted Curve")
				theta[b,i,j] = 1./popt[0]
				omega[b,i,j] = popt[1]
				x = sp.symbols('x')
				print('solving')
				f = interp1d(func_y(t,*popt),t)
				theta[b,i,j] = f(np.exp(-1))
				#print(theta[b,i,j])
				#exit()
				#it = 0
				#while (func_y(it,*popt)>np.exp(-1)):
				#	it +=1
				#theta[b,i,j] = it
				#theta[b,i,j] = sp.solve([sp.exp(-popt[0]*x)*sp.cos(popt[1]*x)-sp.exp(-1),x>0,x<100],x)
				#print(theta[b,i,j])
				#for tau in t:
				#	if func_y(tau,*popt) <= np.exp(-1):
				#		theta[b,i,j] = tau
				#		break
				

			ax[k].plot(t,R[2,b,i,j,:],'b-',label = "Original FFE Curve")
			
			ax[k].grid()
			ax[k].text(.5,.95,r'$ T_L =$ %.2f' % theta[b,i,j], horizontalalignment='center',verticalalignment='center', transform=ax[k].transAxes)
			
			if i == 0:
				ax[k].set_ylim([-.2,1])
				
				
			else:
				ax[k].set_ylim([-.75,1])

			k+=1
	if (b == 0 or b == 6):
		ax[0].legend(loc = 'lower left')
	
	ax[1].set_xlabel('Time Lag (days)')
	ax[3].set_xlabel('Time Lag (days)')
	ax[0].set_ylabel('Top Layer')
	ax[1].set_ylabel('Bottom Layer') 
	ax[0].set_title('R$_x$')
	ax[2].set_title('R$_y$')
	fig_name = fig_dir + 'fitted_pseudo_osc_R_bin%i' % (b+1)
	plt.savefig(fig_name)
	plt.close()
	

	
# WRITE THETA TO FILE

#file_name = home_dir + '/STATS/THETA/pseudo_theta_osc.nc'
#theta_data = Dataset(file_name,'w',format='NETCDF4_CLASSIC')
#theta_data.createDimension('Bin',nbins)
#theta_data.createDimension('Dimension',2)
#theta_data.createDimension('Layer',2)
#theta_var = theta_data.createVariable('Theta',np.float64,('Bin','Dimension','Layer',))
#theta_var[:] = theta 

#file_name = home_dir + '/STATS/THETA/omega.nc'
#omega_data = Dataset(file_name,'w',format = 'NETCDF4_CLASSIC')
#omega_data.createDimension('Bin',nbins)
#omega_data.createDimension('Dimension',2)
#omega_data.createDimension('Layer',2)
#omega_var = omega_data.createVariable('Theta',np.float64,('Bin','Dimension','Layer'))
#omega_var[:] = omega





	
	




