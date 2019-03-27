#!/usr/bin/env python

## CODE FINDS PROPERTIES OF PROPAGATAING EOF MODES ##

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import subprocess
from matplotlib.animation import FuncAnimation

import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')

from calculate_var import velocity_from_psi

regime = 1
ii = 512
jj = 512

# Read the zonal eof pattern from file

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/1/EOF/'

eof_file = home_dir + '/STATS/EOF/eof.nc'
eof_data = Dataset(eof_file,'r')
eof = np.transpose(eof_data.variables['EOFs'][:,:,:,:5])
print(eof.shape)

pc_file = home_dir + '/STATS/EOF/pc.nc'
pc_data = Dataset(pc_file,'r')
pc = np.transpose(pc_data.variables['PCs'][:5,:])
print(pc.shape)
nt = 1000
modes = [4]
psi = np.zeros((ii,jj,2,nt))

for m in modes:
	print(m)
	for t in range(nt):
		for l in range(2):
			psi[:,:,l,t] = psi[:,:,l,t] + eof[m,:,:,l]*pc[t,m]
			
xmin,xmax = 0,520
xx,yy = np.mgrid[xmin:xmax:512j,xmin:xmax:512j]


u,v = velocity_from_psi(psi[:,:,:,-1])

fig = plt.figure()

ax1 = plt.subplot(111)

cax1 = ax1.pcolor(xx,yy,psi[:,:,0,0], cmap = cm.jet,vmin = -50,vmax = 50)

fig.colorbar(cax1)

def animate(i):   
	cax1.set_array(psi[:-1,:-1,0,20*i].flatten())
	return cax1,
        
anim = FuncAnimation(fig, animate, interval=100, frames=int(nt/20))

anim.save('PV.mp4',metadata={'artist':'Guido'},writer='ffmpeg',fps=8,bitrate=500)

plt.show()
plt.close()


plt.pcolor(xx,yy,psi[:,:,0,-1],cmap = cm.jet)
plt.colorbar()
plt.show()

plt.close()

plt.pcolor(xx,yy,u[:,:,0],cmap=cm.jet)
plt.colorbar()
plt.show()
plt.close()

plt.plot(pc)
plt.show()
plt.close()


