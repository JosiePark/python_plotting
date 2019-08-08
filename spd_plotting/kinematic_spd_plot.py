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

regime = 2

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/SPD/KINEMATIC/' % regime
nbins = 10
nt = 1000
nrel = 9

# READ KINEMATIC SPD #

SPD = np.zeros((2,nbins,2,nt))
U0_SPD = np.zeros((2,nbins,2,nt))

file_name = home_dir + '/STATS/SPD/KINEMATIC/kinematic_top.nc'
spd_data = Dataset(file_name,'r')
print(spd_data)
SPD[0,:,:,:] = np.transpose(spd_data.variables['SPD'][:nt,:,:])

file_name = home_dir + '/STATS/SPD/KINEMATIC/kinematic_zonalRossby_top.nc'
spd_data = Dataset(file_name,'r')
print(spd_data)
U0_SPD[0,:,:,:] = np.transpose(spd_data.variables['SPD'][:nt,:,:])


file_name = home_dir + '/STATS/SPD/KINEMATIC/kinematic_bottom.nc'
spd_data = Dataset(file_name,'r')
print(spd_data)
SPD[1,:,:,:] = np.transpose(spd_data.variables['SPD'][:nt,:,:])

file_name = home_dir + '/STATS/SPD/KINEMATIC/kinematic_zonalRossby_bottom.nc'
spd_data = Dataset(file_name,'r')
print(spd_data)
U0_SPD[1,:,:,:] = np.transpose(spd_data.variables['SPD'][:nt,:,:])



# READ PSEUDO SPD

pseudo_SPD = np.zeros((2,nbins,2,nt))

for b in range(nbins):
	file_name = home_dir + '/STATS/SPD/PSEUDO/pseudo_uniform_SPD_bin%i.nc' % (b+1)
	spd_data = Dataset(file_name,'r')
	print(np.transpose(spd_data.variables['SPD'][:,0,:,:]).shape)
	pseudo_SPD[0,b,:,:] = np.mean(np.transpose(spd_data.variables['SPD'][:,0,:,:]),1)
	pseudo_SPD[1,b,:,:] = np.mean(np.transpose(spd_data.variables['SPD'][:,1,:,:]),1)

# PLOT SPD

t = np.arange(0,nt,1)

right_space = 0.02
left_space = 0.08
top_space = 0.05
bottom_space = 0.08
hor_space = 0.06
ver_space = 0.06
nrows = 2
ncols = 2
fig_width = 8.27

for b in range(nbins):
	ax = spd_grid_plot(left_space,right_space,bottom_space,top_space,ver_space,hor_space)
	k=0
	for i in range(ncols):
		for j in range(nrows):
			ax[k].plot(t,SPD[j,b,i,:],'b-',label = 'rossbyHalf')
			ax[k].plot(t,U0_SPD[j,b,i,:],'g-.',label = 'rossbyHalf +  zonal')
			ax[k].plot(t,pseudo_SPD[j,b,i,:],'r--',label = 'FFE')
			if (i==0 and j == 0):
				ax[k].legend(loc = 'upper left')
			ax[k].ticklabel_format(style = 'sci',axis = 'y',scilimits=(0,0))
			ax[k].grid()
			
			if j == 0 and i == 0:
				ax[k].set_title('D$_x$ (km s $^{-1}$)')
			elif j == 0 and i == 1:
				ax[k].set_title('D$_y$ (km s $^{-1}$)')
			
			if j == 0 and i == 0:
				ax[k].set_ylabel('Top Layer')
			elif j == 1 and i == 0:
				ax[k].set_ylabel('Bottom Layer')
				
			if j == 1:
				ax[k].set_xlabel('Time (Days)')
				
			if j == 0:
				if i == 0:
					ax[k].set_ylim([0.,8.e6])
				else:
					ax[k].set_ylim([0.,2.e4])
			else:
				if i == 0:
					ax[k].set_ylim([0.,8.e5])
				else:
					ax[k].set_ylim([0.,2.5e4])
			k+=1	
		
	#ax[0].set_ylim([0,3.e4])
	#ax[1].set_ylim([0,5])	

	fig_name = fig_dir + 'zonal_SPD_bin%i' % (b+1)
	plt.savefig(fig_name)
	
	plt.close()

# compare contributions across the bins

# READ TIME SCALE #

file_name = home_dir + '/STATS/TSCALE/pseudo_MEAN_TSCALE.nc'
tscale_data = Dataset(file_name,'r')
tscale = np.transpose(tscale_data.variables['Time Scale'][:])

ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

contribution = np.zeros((2,2,nbins,2))

for l in range(2):
	for b in range(nbins):
		if int(tscale[b,l]) == 1000:
			T = nt-1
		else:
			T = int(tscale[b,l])-1
			
		contribution[0,:,b,l] = (pseudo_SPD[:,b,l,T]-SPD[:,b,l,T])/pseudo_SPD[:,b,l,T]
		contribution[1,:,b,l] = (pseudo_SPD[:,b,l,T]-U0_SPD[:,b,l,T])/pseudo_SPD[:,b,l,T]


k = 0

for c in range(ncols):
	for r in range(nrows):
		line = []
		ax[k].plot(contribution[0,c,:,r],np.linspace(1,nbins,nbins),label = 'rossbyHalf')
		ax[k].plot(contribution[1,c,:,r],np.linspace(1,nbins,nbins),label = 'rossbyHalf + zonal')
		
		ax[k].set_xlim([-2,2])
		if (r == 1):
			ax[k].set_xlabel('Contribution to SPD')
			if (c == 0):
				ax[k].set_ylabel('Bottom Layer')
		if (r == 0) and (c == 0):
			ax[k].set_ylabel('Top Layer')
			ax[k].set_title('Zonal')
			ax[k].legend(loc = 'upper left')
		if (r == 0) and (c == 1):
			ax[k].set_title('Meridional')
		if (c==1):
			ax[k].set_yticklabels([])
		if (r == 0):
			ax[k].set_xticklabels([])
		ax[k].grid()
		
	

		k+=1

fig_name = fig_dir + 'contribution'
plt.savefig(fig_name)

plt.close()

#plt.show()

################################################################

## Examine effect of altering the amplitude

## READ DATA

SPD = np.zeros((3,2,nbins,2,nt))

file_name = home_dir + '/STATS/SPD/KINEMATIC/kinematic_top.nc'
spd_data = Dataset(file_name,'r')
SPD[0,0,:,:,:] = np.transpose(spd_data.variables['SPD'][:nt,:,:])

file_name = home_dir + '/STATS/SPD/KINEMATIC/kinematic_bottom.nc'
spd_data = Dataset(file_name,'r')
SPD[0,1,:,:,:] = np.transpose(spd_data.variables['SPD'][:nt,:,:])

file_name = home_dir + '/STATS/SPD/KINEMATIC/kinematic_top_2A.nc'
spd_data = Dataset(file_name,'r')
SPD[1,0,:,:,:] = np.transpose(spd_data.variables['SPD'][:nt,:,:])

file_name = home_dir + '/STATS/SPD/KINEMATIC/kinematic_bottom_2A.nc'
spd_data = Dataset(file_name,'r')
SPD[1,1,:,:,:] = np.transpose(spd_data.variables['SPD'][:nt,:,:])

file_name = home_dir + '/STATS/SPD/KINEMATIC/kinematic_top_halfA.nc'
spd_data = Dataset(file_name,'r')
SPD[2,0,:,:,:] = np.transpose(spd_data.variables['SPD'][:nt,:,:])

file_name = home_dir + '/STATS/SPD/KINEMATIC/kinematic_bottom_halfA.nc'
spd_data = Dataset(file_name,'r')
SPD[2,1,:,:,:] = np.transpose(spd_data.variables['SPD'][:nt,:,:])

for b in range(nbins):
	ax = spd_grid_plot(left_space,right_space,bottom_space,top_space,ver_space,hor_space)
	k=0
	for i in range(ncols):
		for j in range(nrows):
			ax[k].plot(t,SPD[0,j,b,i,:],'b-',label = 'A')
			ax[k].plot(t,SPD[1,j,b,i,:],'g-.',label = '2A')
			ax[k].plot(t,SPD[2,j,b,i,:],'m-',label = 'A/2')
			ax[k].plot(t,pseudo_SPD[j,b,i,:],'r--',label = 'FFE')
			if (i==0 and j == 0):
				ax[k].legend(loc = 'upper left')
			ax[k].ticklabel_format(style = 'sci',axis = 'y',scilimits=(0,0))
			ax[k].grid()
			
			if j == 0 and i == 0:
				ax[k].set_title('D$_x$ (km s $^{-1}$)')
			elif j == 0 and i == 1:
				ax[k].set_title('D$_y$ (km s $^{-1}$)')
			
			if j == 0 and i == 0:
				ax[k].set_ylabel('Top Layer')
			elif j == 1 and i == 0:
				ax[k].set_ylabel('Bottom Layer')
				
			if j == 1:
				ax[k].set_xlabel('Time (Days)')
				
			if j == 0:
				if i == 0:
					ax[k].set_ylim([0.,8.e6])
				else:
					ax[k].set_ylim([0.,2.e4])
			else:
				if i == 0:
					ax[k].set_ylim([0.,8.e5])
				else:
					ax[k].set_ylim([0.,2.5e4])
			k+=1	
		
	#ax[0].set_ylim([0,3.e4])
	#ax[1].set_ylim([0,5])	

	fig_name = fig_dir + 'amplitude_SPD_bin%i' % (b+1)
	plt.savefig(fig_name)
	
	plt.close()




	
	




