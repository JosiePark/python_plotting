## CODE THAT CALCULATES PROPAGAION SPEED OF EOFS 1-2

#!/usr/bin/env python

## CODE FINDS PROPERTIES OF PROPAGATAING EOF MODES ##

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')

from calculate_var import velocity_from_psi
from grid_plot import square_grid_plot

# READ FIRST TWO EOFS

regime = 1
ii = 512
jj = 512
nt = 400
	
home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/1/EOF/'

field_file = home_dir + '/QG/eof_1-2.nc'
field_data = Dataset(field_file,'r')
psi = np.transpose(field_data.variables['Stream Function'][:nt,0,:,:])

eof_file = home_dir + '/STATS/EOF/eof.nc'
eof_data = Dataset(eof_file,'r')
eof = np.transpose(eof_data.variables['EOFs'][0,:,:,:2])

pc_file = home_dir + '/STATS/EOF/pc.nc'
pc_data = Dataset(pc_file,'r')
pc = np.transpose(pc_data.variables['PCs'][:,:2])

## out of phase by a day

## work out spatial phase difference between two eofs

## plot hovmoler diagram


xcoord = 200
eps = 10
max_ycoord = np.argmax(psi[xcoord,:,10])
center_start = np.argmax(psi[xcoord,max_ycoord,:40])
center_end = np.argmax(psi[xcoord,max_ycoord,center_start+eps:])
print(max_ycoord,center_start,center_end+eps+center_start)

plt.pcolor(psi[xcoord,:,:],cmap = cm.jet)
plt.xlabel('Time (days)')
plt.ylabel('Y (km)')
plt.ylim([0,520])
plt.colorbar()
fig_name = fig_dir + 'hovmoller_eof_1-2'
plt.savefig(fig_name)







