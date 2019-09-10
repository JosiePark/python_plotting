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
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/SPD/STOCHASTIC_MODELS/DIFFUSION/' % regime
nbins = 10
nt = 1000

# READ SPD #

SPD = np.zeros((3,nbins,nt))

for b in range(nbins):
	file_name = home_dir + '/STATS/SPD/FULL/full_uniform_SPD_bin%i.nc' % (b+1)
	spd_data = Dataset(file_name,'r')
	tmp = np.transpose(spd_data.variables['SPD'][:])
	#print(tmp.shape)
	#print(spd_data)
	SPD[0,b,:] = np.mean(tmp[1,:,0,:nt],0)
	
	file_name = home_dir + '/STATS/SPD/FINAL_DIFFUSION/PV_BIN/pv_bin_diffusion_SPD_bin%i.nc' % (b+1)
	spd_data = Dataset(file_name,'r')
	tmp = np.transpose(spd_data.variables['SPD'][:,0,1])
	#print(spd_data)
	SPD[1,b,:] = tmp[:nt]

	file_name = home_dir + '/STATS/SPD/FINAL_DIFFUSION/PV_BIN_PV/pv_bin_pv_diffusion_SPD_bin%i.nc' % (b+1)
	spd_data = Dataset(file_name,'r')
	tmp = np.transpose(spd_data.variables['SPD'][:])
	#print(spd_data)
	SPD[2,b,:] = tmp[:nt]
	
	#exit()
			
# PLOT SPD

t = np.arange(0,1000,1)
print(SPD[2,0,2])

hor_space = 0.0
ver_space = 0.02
nrows = 5
ncols = 2
left_space = 0.1
right_space = 0.01
bottom_space = 0.05
top_space = 0.02

ax = a4_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,ver_space,hor_space)



k = 0

#SPD[2,:,1] = 0


for i in range(ncols):
	for j in range(nrows):
		ax[k].plot(t,SPD[0,k,:],'b-',label = 'Full')
		ax[k].plot(t,SPD[1,k,:],'r--',label = 'Diffusion')
		ax[k].plot(t,SPD[2,k,:],'g-.',label = 'PV Mapped Diffusion')
		
		ax[k].ticklabel_format(style = 'sci',axis = 'y',scilimits=(0,0))	
		#ax[k].set_ylim([0,4.e7])
		#ax[k].set_ylim([0,1.5e8])
		ax[k].set_ylim([0,1.e4])
		ax[k].text(.9,.1,'Bin %i' % (k+1),horizontalalignment = 'center',verticalalignment = 'center', transform = ax[k].transAxes)
		if i == 0:
			ax[k].set_ylabel('D$_y$ (km $^{2}$)')
		else:
			ax[k].set_yticklabels([])
		if j == nrows-1:
			ax[k].set_xlabel('Time (days)')
		else:
			ax[k].set_xticklabels([])
		ax[k].grid()
		k+=1
		

ax[0].legend()
fig_name = fig_dir + 'pvbin_pv_diffusion_SPD'
plt.savefig(fig_name)
plt.close()
ax = a4_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,ver_space,hor_space)
k = 0

for i in range(ncols):
	for j in range(nrows):
		ax[k].loglog(t,SPD[0,k,:],'b-',label = 'Full')
		ax[k].loglog(t,SPD[1,k,:],'r--',label = 'Diffusion')
		ax[k].loglog(t,SPD[2,k,:],'g-.',label = 'PV Mapped Diffusion')
		ax[k].loglog(t[:100],t[:100]**2,'k--',label = 'Ballisitc')
		ax[k].loglog(t[100:],t[100:],'k-.',label = 'Diffusive')
		
		#ax[k].ticklabel_format(style = 'sci',axis = 'y',scilimits=(0,0))	
		#ax[k].set_ylim([1,1.e5])
		ax[k].text(.9,.1,'Bin %i' % (k+1),horizontalalignment = 'center',verticalalignment = 'center', transform = ax[k].transAxes)
		if i == 0:
			ax[k].set_ylabel('D$_y$ (km $^{2}$)')
		else:
			ax[k].set_yticklabels([])
		if j == nrows-1:
			ax[k].set_xlabel('Time (days)')
		else:
			ax[k].set_xticklabels([])
		ax[k].grid()
		k+=1
		

ax[0].legend()
fig_name = fig_dir + 'pvbin_pv_diffusion_logSPD'
plt.savefig(fig_name)
	
	




