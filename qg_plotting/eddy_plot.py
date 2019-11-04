#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import sys

sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import square_grid_plot
from calculate_var import pv_anomaly_from_psi,velocity_from_psi

sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/input/')
from parameters import *

# INPUT PARAMETERS #

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/'
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/COMBO/QG/'


# READ PSI #

qg_file = home_dir + '1/QG/QG.nc'
qg_data = Dataset(qg_file,'r')
tmp = np.transpose(qg_data.variables['Stream Function'][0,:,:,:])
[ii,jj,nlayers] = tmp.shape
psi = np.zeros((2,ii,jj,nlayers))
psi[0,:,:,:] = tmp

qg_file = home_dir + '2/QG/QG.nc'
qg_data = Dataset(qg_file,'r')
tmp = np.transpose(qg_data.variables['Stream Function'][0,:,:,:])
psi[1,:,:,:] = tmp

# READ TIME-AVERAGED STREAM FUNCTCION

psi_av = np.zeros((2,ii,jj,nlayers))

ave_file = home_dir + '1/QG/QG_ave.nc'
ave_data = Dataset(ave_file,'r')
psi_av[0,:,:,:] = np.transpose(ave_data.variables['Time-averaged Stream Function'][:,:,:])


ave_file = home_dir + '2/QG/QG_ave.nc'
ave_data = Dataset(ave_file,'r')
psi_av[1,:,:,:] = np.transpose(ave_data.variables['Time-averaged Stream Function'][:,:,:])

psi = psi-psi_av


# CALCLATE PV ANOMALY

zeta = np.zeros((psi.shape))
zeta[0,:,:,:] = pv_anomaly_from_psi(psi[0,:,:,:],H1,H2,Rd,basinscale)
zeta[1,:,:,:] = pv_anomaly_from_psi(psi[1,:,:,:],H1,H2,Rd,basinscale)

# PLOT PV ANOMALY #

left_space = .08
top_space = .05
bottom_space = .08
right_space = .05
nrows = 2
ncols = 2
fig_width = 8.27
hor_space = .04
ver_space = .02

ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

xmin,xmax = 0,520
xx,yy = np.mgrid[xmin:xmax:512j,xmin:xmax:512j]

k = 0

for i in range(ncols):
	for j in range(nrows):
		c = ax[k].pcolor(xx,yy,psi[i,:,:,j],cmap = cm.jet)
		if j == 1:
			ax[k].set_xlabel('X (km)')
		else:
			ax[k].set_xticklabels([])
		if i == 0:
			ax[k].set_ylabel('Y (km)')
		else:
			ax[k].set_yticklabels([])

		plt.colorbar(c,ax=ax[k])
		

		k+=1
		
ax[0].set_title('Coherent Jet')
ax[2].set_title('Latent Jet')
fig_name = fig_dir + 'eddy_psi'
plt.savefig(fig_name)
plt.show()

plt.close()

## CALCULATE VELOCTY OF eddying STREAM FUNCTION

u_eddy = np.zeros((2,2,ii,jj,2))

for i in range(2):

	[u_eddy[i,0,:,:,:],u_eddy[i,1,:,:,:]] = velocity_from_psi(psi[i,:,:,:])
	
ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

k = 0

for i in range(ncols):
	for j in range(nrows):
		c = ax[k].pcolor(xx,yy,u_eddy[i,0,:,:,j],cmap = cm.jet)
		if j == 1:
			ax[k].set_xlabel('X (km)')
		else:
			ax[k].set_xticklabels([])
		if i == 0:
			ax[k].set_ylabel('Y (km)')
		else:
			ax[k].set_yticklabels([])

		cbar = plt.colorbar(c,ax=ax[k])
		cbar.ax.set_ylabel('u\' (cm s$^{-1}$)')

		k+=1
		
ax[0].set_title('Coherent Jet')
ax[2].set_title('Latent Jet')
fig_name = fig_dir + 'eddy_u'
plt.savefig(fig_name)
plt.show()
plt.close()

ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

k = 0

for i in range(ncols):
	for j in range(nrows):
		c = ax[k].pcolor(xx,yy,u_eddy[i,1,:,:,j],cmap = cm.jet)
		if j == 1:
			ax[k].set_xlabel('X (km)')
		else:
			ax[k].set_xticklabels([])
		if i == 0:
			ax[k].set_ylabel('Y (km)')
		else:
			ax[k].set_yticklabels([])

		cbar = plt.colorbar(c,ax=ax[k])
		cbar.ax.set_ylabel('v\' (cm s$^{-1}$)')

		k+=1
		
ax[0].set_title('Coherent Jet')
ax[2].set_title('Latent Jet')
fig_name = fig_dir + 'eddy_v'
plt.savefig(fig_name)
plt.show()
plt.close()

ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

k = 0

for i in range(ncols):
	for j in range(nrows):
		c = ax[k].pcolor(xx,yy,psi_av[i,:,:,j],cmap = cm.jet)
		if j == 1:
			ax[k].set_xlabel('X (km)')
		else:
			ax[k].set_xticklabels([])
		if i == 0:
			ax[k].set_ylabel('Y (km)')
		else:
			ax[k].set_yticklabels([])

		plt.colorbar(c,ax=ax[k])
		ax[0].set_title('Coherent Jet')
		ax[2].set_title('Latent Jet')

		k+=1
fig_name = fig_dir + 'mean_psi'
plt.savefig(fig_name)
plt.show()
plt.close()

## CALCULATE VELOCTY OF TIME AVERAGED STREAM FUNCTION

u_av = np.zeros((2,2,ii,jj,2))

for i in range(2):

	[u_av[i,0,:,:,:],u_av[i,1,:,:,:]] = velocity_from_psi(psi_av[i,:,:,:])

ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

k = 0

for i in range(ncols):
	for j in range(nrows):
		c = ax[k].pcolor(xx,yy,u_av[i,0,:,:,j],cmap = cm.jet)
		if j == 1:
			ax[k].set_xlabel('X (km)')
		else:
			ax[k].set_xticklabels([])
		if i == 0:
			ax[k].set_ylabel('Y (km)')
		else:
			ax[k].set_yticklabels([])

		cbar = plt.colorbar(c,ax=ax[k])
		cbar.ax.set_ylabel('U (cm s$^{-1}$)')

		k+=1
		
ax[0].set_title('Coherent Jet')
ax[2].set_title('Latent Jet')
fig_name = fig_dir + 'mean_u'
plt.savefig(fig_name)
plt.show()
plt.close()




	
	
	




