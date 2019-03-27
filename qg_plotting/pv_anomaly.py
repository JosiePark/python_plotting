#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import sys

sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import square_grid_plot
from calculate_var import pv_anomaly_from_psi, velocity_from_psi

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
		c = ax[k].pcolor(xx,yy,zeta[i,:,:,j],cmap = cm.jet)
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
fig_name = fig_dir + 'PV_anomaly'
plt.savefig(fig_name)
plt.show()

## PLOT ZONALLY AVERAGED VELOCITY

u = np.zeros((2,ii,jj,2))
v = np.zeros((2,ii,jj,2))
for i in range(2):

	[u[i,:,:,:],v[i,:,:,:]] = velocity_from_psi(psi[i,:,:,:])

u_zonal = np.mean(u,axis = 1)

plt.close()

y = np.arange(0.,512.,1)

k = 0

hor_space = .02
ver_space = .05

ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

for i in range(ncols):
	for j in range(nrows):
		c = ax[k].plot(u_zonal[i,:,j],y)
		if j == 1:
			ax[k].set_xlabel('Zonal Velocity')
		if i == 0:
			ax[k].set_ylabel('Y (km)')
		else:
			ax[k].set_yticklabels([])
		ax[0].set_title('Coherent Jet')
		ax[2].set_title('Latent Jet')

		k+=1
fig_name = fig_dir + 'zonal_velocity'
plt.savefig(fig_name)
plt.show()
	







	
	
	




