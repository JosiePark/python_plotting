#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import spd_grid_plot

# INPUT PARAMETERS #

regime = 2

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/SPD/STOCHASTIC_MODELS/MARKOV1/' % regime
nbins = 10
nt = 1000

# READ SPD #

SPD = np.zeros((3,nbins,2,2,nt))

for b in range(nbins):
	file_name = home_dir + '/STATS/SPD/FULL/full_uniform_SPD_bin%i.nc' % (b+1)
	spd_data = Dataset(file_name,'r')
	tmp = np.transpose(spd_data.variables['SPD'][:])
	SPD[0,b,:,:,:] = np.mean(tmp[:,:,:,:nt],1)
	
	file_name = home_dir + '/STATS/SPD/FINAL_MARKOV1/DIFF/diff_markov1_SPD_bin%i.nc' % (b+1)
	spd_data = Dataset(file_name,'r')
	tmp = np.transpose(spd_data.variables['SPD'][:])
	SPD[1,b,:,:,:] = tmp[:,:,:nt]

	file_name = home_dir + '/STATS/SPD/FINAL_MARKOV1/LOOPING/looping_markov1_SPD_bin%i.nc' % (b+1)
	spd_data = Dataset(file_name,'r')
	tmp = np.transpose(spd_data.variables['SPD'][:])
	SPD[2,b,:,:,:] = tmp[:,:,:nt]
	
	#exit()
			
# PLOT SPD

t = np.arange(0,1000,1)

nrows = 2
ncols = 2
right_space = 0.02
left_space = 0.12
top_space = 0.05
bottom_space = 0.09
hor_space = 0.08
ver_space = 0.08
fig_width = 8.27

for b in range(nbins):
	ax = spd_grid_plot(left_space,right_space,bottom_space,top_space,ver_space,hor_space)
	k=0
	for i in range(ncols):
		for j in range(nrows):
			ax[k].plot(t,SPD[0,b,i,j,:],'b-',label = 'Full')
			ax[k].plot(t,SPD[1,b,i,j,:],'g--',label = 'Markov1')
			ax[k].plot(t,SPD[2,b,i,j,:],'r-.',label = 'Looping Markov1')
			#ax[k].axvline(x = tscale[b,j],color = 'k',label = 'T')
			#ax[k].axvline(x = theta[b,i,j], color = 'k', linestyle = '-.',label = 'T_L')

			ax[k].ticklabel_format(style = 'sci',axis = 'y',scilimits=(0,0))
			ax[k].grid()
			k+=1
	if (regime == 1):
		ax[0].set_ylim([0,8.e7])
		ax[1].set_ylim([0,2.e7])
		ax[2].set_ylim([0,2.5e4])
		ax[3].set_ylim([0,1.2e4])
		
	else:
		ax[0].set_ylim([0,1.e8])
		ax[1].set_ylim([0,1.e7])
		ax[2].set_ylim([0,2.e4])
		ax[3].set_ylim([0,2.e4])
	
	if (b == 0 or b == 6):
		ax[0].legend(loc = 'upper left')
	ax[1].set_xlabel('Time (days)')
	ax[3].set_xlabel('Time (days)')
	ax[0].set_ylabel('Top Layer')
	ax[1].set_ylabel('Bottom Layer') 
	ax[0].set_title('D$_x$ (km$^{2}$)')
	ax[2].set_title('D$_y$ (km$^{2}$)')
	fig_name = fig_dir + 'new_looping_markov1_SPD_bin%i' % (b+1)
	plt.savefig(fig_name)
	
for b in range(nbins):
	ax = spd_grid_plot(left_space,right_space,bottom_space,top_space,ver_space,hor_space)
	k=0
	for i in range(ncols):
		for j in range(nrows):
			ax[k].loglog(t,SPD[0,b,i,j,:],'b-',label = 'Full')
			ax[k].loglog(t,SPD[1,b,i,j,:],'g--',label = 'Markov1')
			ax[k].loglog(t,SPD[2,b,i,j,:],'r-.',label = 'Looping Markov1')
			ax[k].loglog(t[:100],t[:100]**2,'k--',label = 'Ballisitc')
			ax[k].loglog(t[100:],t[100:],'k-.',label = 'Diffusive')
			#ax[k].axvline(x = tscale[b,j],color = 'k',label = 'T')
			ax[k].axvline(x = theta[b,i,j], color = 'k', linestyle = '-',label = '$T_L$')

			#ax[k].ticklabel_format(style = 'sci',axis = 'y',scilimits=(0,0))
			ax[k].grid()
			k+=1
	if (b == 0 or b == 6):
		ax[0].legend(loc = 'lower left')
	#if (regime == 1):
	#	ax[0].set_ylim([0,8.e7])
	#	ax[1].set_ylim([0,2.e7])
	#	ax[2].set_ylim([0,2.5e4])
	#	ax[3].set_ylim([0,1.2e4])
		
	#else:
	#	ax[0].set_ylim([0,1.e8])
	#	ax[1].set_ylim([0,1.e7])
	#	ax[2].set_ylim([0,1.e5])
	#	ax[3].set_ylim([0,2.e4])
	
			
	ax[1].set_xlabel('Time (days)')
	ax[3].set_xlabel('Time (days)')
	ax[0].set_ylabel('Top Layer')
	ax[1].set_ylabel('Bottom Layer') 
	ax[0].set_title('D$_x$ (km$^{2}$)')
	ax[2].set_title('D$_y$ (km$^{2}$)')
	fig_name = fig_dir + 'new_looping_markov1_logSPD_bin%i' % (b+1)
	plt.savefig(fig_name)
	




