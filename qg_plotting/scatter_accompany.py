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

regime = 2
home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i/' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/QG/' % regime


# READ PSI #

read_time = [0,50,100,150]
qg_file = home_dir + 'QG/QG.nc'
qg_data = Dataset(qg_file,'r')
psi = np.transpose(qg_data.variables['Stream Function'][read_time,:,:,:])
time = qg_data.variables['Time'][:]
[ii,jj,nlayers,nt] = psi.shape
print(time[read_time])


# CALCLATE PV ANOMALY

zeta = np.zeros((psi.shape))
for t in range(nt):
	zeta[:,:,:,t] = pv_anomaly_from_psi(psi[:,:,:,t],H1,H2,Rd,basinscale)

# PLOT PV ANOMALY #

left_space = .05
top_space = .05
bottom_space = .08
right_space = .04
nrows = 2
ncols = 2
fig_width = 8.27
hor_space = .06
ver_space = .08

ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

xmin,xmax = 0,520
xx,yy = np.mgrid[xmin:xmax:512j,xmin:xmax:512j]

read = [0, 2, 1, 3]
for k in range(4):
	c = ax[k].pcolor(xx,yy,zeta[:,:,0,read[k]],cmap = cm.jet)
	ax[k].set_title('Time = %i days' % round(time[read_time[read[k]]]))
	plt.colorbar(c,ax=ax[k])

fig_name = fig_dir + 'PV_anomaly_snapshot'
plt.savefig(fig_name)

	







	
	
	




