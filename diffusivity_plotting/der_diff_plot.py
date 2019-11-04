#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import square_grid_plot,a4_plot,spd_grid_plot

plt.rcParams.update({'font.size': 8})

# INPUT PARAMETERS #

regime = 1

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i/' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/STATS/DIFFUSIVITY/' % regime
nbins = 10
nt = 1000

# READ TIME-AVERAGED STREAM FUNCTION

file_name = home_dir + 'QG/QG_ave.nc'
psi_data = Dataset(file_name,'r')
psi_ave = np.transpose(psi_data.variables['Time-averaged Stream Function'][:])

# READ DIFFUSIVITY

K = np.zeros((2,nbins,2))
	
file_name = home_dir + 'STATS/DIFFUSIVITY/pseudo_new_DIFF.nc'
diff_data = Dataset(file_name,'r')
tmp = np.transpose(diff_data.variables['Diffusivity'][:])
print(diff_data)
K = tmp

# READ LAGRANGIAN VELOCITY VARIANCE

file_name = home_dir + '/STATS/SIGMA/new_lagrangian_sigma.nc'
sigma_data = Dataset(file_name,'r')
sigma_L = sigma_data.variables['Lagrangian Velocity Variance'][:]
print(sigma_data)
plt.plot(sigma_L[:,0,0,0])
plt.show()
sigma = np.zeros((2,10,2))
sigma[0,:,:] = sigma_L[:,0,0,:]
sigma[1,:,:] = sigma_L[:,1,1,:]

theta_dir = home_dir + 'STATS/THETA/'
theta_file = theta_dir + 'pseudo_theta_osc.nc'
theta_data = Dataset(theta_file,'r')
print(theta_data)
theta = theta_data.variables['Theta'][:]
theta = np.swapaxes(theta,0,1)

# RESCALE SIGMA	

basinscale = 520.e5
ii = 512
scale = basinscale/float(ii)
uscale = 1.
tscale = scale/uscale

scale_new = scale/(10**5)
tscale_new = tscale/86400.

#sigma = sigma*scale_new**2

K_new = np.multiply(sigma,theta)

K_new = K_new/86400.

# PLOT ALPHA SUPERIMPOSED ON THE TIME-AVERAGED STREAM FUNCTION

xmin,xmax = 0,520

xx,yy = np.mgrid[xmin:xmax:512j,xmin:xmax:512j]

bin_width = 512./10.
bin_centres = []
bin_centres.append(bin_width/2.)
for b in range(nbins-1):
	bin_centres.append(bin_centres[b]+bin_width)
	
right_space = 0.02
left_space = 0.1
top_space = 0.08
bottom_space = 0.08
hor_space = 0.07
ver_space = 0.06
nrows = 2
ncols = 2
fig_width = 8.27

k = 0
ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

K = K/86400.

K = K*(1000**2)
K_new = K_new*(1000**2)

for i in range(ncols):
	for j in range(nrows):
		
		ax[k].pcolor(xx,yy,psi_ave[:,:,j],alpha=0.5,cmap = cm.gray)
		ax_K = ax[k].twiny()
		ax_K.plot(K[i,:,j],bin_centres,'b-',label='From D')
		ax_K.plot(K_new[i,:,j],bin_centres,'g--',label='From sigma*T')
		if (j == 0):
			if (i == 0):
				ax_K.set_xlabel('K$_x$ (m$^{2}$ s$^{-1}$)')
			else:
				ax_K.set_xlabel('K$_y$ (m$^{2}$ s$^{-1}$)')
			ax[k].set_xticklabels([])
		else:
			ax[k].set_xlabel('X (km)')
		#ax_K.set_xlim(0,3)
		ax[k].set_ylim(0,520)
		ax[k].set_xlim(0,520)
		#ax_K.ticklabel_format(style = 'sci',axis = 'x',scilimits=(0,0))
		if (i==0):
			ax[k].set_ylabel('Y (km)')
		else:
			ax[k].set_yticklabels([])
		
		if (i == 0 and j == 0):
			ax_K.legend(loc = 'upper left')
		ax_K.grid()
		k+=1

fig_name = fig_dir + 'test_K_compare'
plt.savefig(fig_name)

	



