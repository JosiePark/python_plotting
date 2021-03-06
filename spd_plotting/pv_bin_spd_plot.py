#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import square_grid_plot,a4_plot

plt.rcParams.update({'font.size': 8})


# INPUT PARAMETERS #

regime = 1

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/SPD/PV_BINS/' % regime
nbins = 10
nrel = 9
nt = 1000

# READ SPD #

SPD = np.zeros((2,nbins,nt))

for b in range(nbins):
	file_name = home_dir + '/STATS/SPD/PV_BINS/TEST_FULL/test_full_PVBINS_SPD_bin%i.nc' % (b+1)
	spd_data = Dataset(file_name,'r')
	tmp = np.transpose(spd_data.variables['SPD'][:])
	SPD[0,b,:] = np.mean(tmp[1,:,0,:],0)
	
file_name = home_dir + '/STATS/PVDISP/test_full_PVDISP.nc' 
spd_data = Dataset(file_name,'r')
tmp = np.transpose(spd_data.variables['PV Mapped Dispersion'][:])
SPD[1,:,:] = np.mean(tmp,1)

# READ TIME SCALE #

file_name = home_dir + '/STATS/TSCALE/full_MEAN_PVBINS_TSCALE.nc'
tscale_data = Dataset(file_name,'r')
tscale = np.transpose(tscale_data.variables['Time Scale'][:])
print(tscale_data)

file_name = home_dir + '/STATS/TSCALE/test_full_PVDISP_MEAN_TSCALE.nc'
tscale_data = Dataset(file_name,'r')
pv_tscale = np.transpose(tscale_data.variables['Time Scale'][:])
print(tscale_data)

print(tscale.shape)
print(pv_tscale.shape)

file_name = home_dir + '/STATS/THETA/theta_PV.nc'
theta_data = Dataset(file_name,'r')
theta = theta_data.variables['Theta'][:]
			
# PLOT SPD

t = np.arange(0,1000,1)

for b in range(nbins):
	fig = plt.figure(constrained_layout = True)
	gs = gridspec.GridSpec(ncols=1,nrows=1,figure=fig)
	ax = []
	ax.append(plt.subplot(gs[0]))
	ax[0].plot(t,SPD[0,b,:],label = 'Full')
	ax[0].plot(t,SPD[1,b,:],label = 'PV Mapped')
	#ax[0].axvline(x = theta[b], label = '$T_L$')
	#ax[0].axvline(x = pv_tscale[b], label = 'PV Mapped T')
	ax[0].legend()
	ax[0].ticklabel_format(style = 'sci',axis = 'y',scilimits=(0,0))
	ax[0].set_xlabel('Time (days)')
	ax[0].set_ylabel('D$_y$ (km$^{2}$)')
	ax[0].grid()
	
	fig_name = fig_dir + 'PV_SPD_bin%i' % (b+1)
	plt.savefig(fig_name)
	

plt.close(fig)

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

for i in range(ncols):
	for j in range(nrows):
		ax[k].plot(t,SPD[0,k,:],'b-',label='Full')
		ax[k].plot(t,SPD[1,k,:],'r-.',label='PV Mapped')
		#ax[k].axvline(x = t[k,0], label = 'T',color = 'k')
		#ax[k].axvline(x = pv_tscale[k], label = 'PV Mapped T',color = 'k',linestyle = '-.')
		ax[k].grid()
		#ax[k].ticklabel_format(style = 'sci',axis = 'y',scilimits=(0,0))
		if regime == 1:
			ax[k].set_ylim([0,6000])
		else:
			ax[k].set_ylim([0,12000])
		if j == nrows-1:
			ax[k].set_xlabel('Time (days)')
		else:
			ax[k].set_xticklabels([])
		if i == 0:
			ax[k].set_ylabel('D$_y$ (km$^{2}$)')
		else:
			ax[k].set_yticklabels([])

		ax[k].text(.9,.05,'Bin %i' % (k+1),horizontalalignment = 'center',verticalalignment = 'center', transform = ax[k].transAxes)
		k+=1
ax[0].legend(loc = 'upper left')

fig_name = fig_dir + 'test_PV_SPD'
plt.savefig(fig_name)
plt.close(fig)

ax = a4_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space)

k = 0

for i in range(ncols):
	for j in range(nrows):
		ax[k].loglog(t,SPD[0,k,:],'b-',label='Full')
		ax[k].loglog(t,SPD[1,k,:],'r-.',label='PV Mapped')
		ax[k].axvline(x = theta[k], label = '$T_L$',color = 'k')
		#ax[k].axvline(x = pv_tscale[k], label = 'PV Mapped T',color = 'k',linestyle = '-.')
		ax[k].grid()
		#ax[k].ticklabel_format(style = 'sci',axis = 'y',scilimits=(0,0))
		#if regime == 1:
		#	ax[k].set_ylim([0,6000])
		#else:
		#	ax[k].set_ylim([0,12000])
		#if j == nrows-1:
		#	ax[k].set_xlabel('Time (days)')
		#else:
		#	ax[k].set_xticklabels([])
		#if i == 0:
		#	ax[k].set_ylabel('D$_y$ (km s $^{-1}$)')
		#else:
			#ax[k].set_yticklabels([])

		ax[k].text(.9,.05,'Bin %i' % (k+1),horizontalalignment = 'center',verticalalignment = 'center', transform = ax[k].transAxes)
		k+=1
ax[0].legend(loc = 'upper left')

fig_name = fig_dir + 'test_PV_logSPD'
plt.savefig(fig_name)
		


	
	




