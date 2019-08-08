#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation

import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')

from calculate_var import pv_anomaly_from_psi

def animate(i):
	cax.set_array(pv[:-1,:-1,0,i].flatten())
	cax.text(.5,.1,'%i days' % int(time[i]))
	return cax
	
regime = 2
ii = 512
jj = 512
H1 = 1.e5
H2 = 3.e5
Rd = 25.e5
basinscale = 520.e5

# read data 

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/QG/' % regime

dt_read = 50

psi_file = home_dir + '/QG/QG.nc'
psi_data = Dataset(psi_file,'r')
psi = np.transpose(psi_data.variables['Stream Function'][0:-1:dt_read,:,:,:])
time = np.transpose(psi_data.variables['Time'][0:-1:dt_read])
ene = np.transpose(psi_data.variables['Energy'][0:-1:dt_read])


# plot energy 

plt.plot(time,np.transpose(ene))
plt.legend(['Potential','Kinetic L1','Kinetic L2','Total'])
plt.xlabel('Time (days)')
plt.ylabel('Energy')
plt.savefig(fig_dir + 'energy')
plt.show()

# calculate pv anomaly
nt = psi.shape[-1]
pv = np.zeros((psi.shape))
for t in range(nt):
	print(t)
	pv[:,:,:,t] = pv_anomaly_from_psi(psi[:,:,:,t],H1,H2,Rd,basinscale)

# animate movie of pv anomaly

xmin,xmax = 0,520
xx,yy = np.mgrid[xmin:xmax:512j,xmin:xmax:512j]

fig = plt.figure()
ax1 = plt.subplot(111)
cax = ax1.pcolor(xx,yy,pv[:,:,0,0], cmap = cm.jet,vmin = -50,vmax = 50)
fig.colorbar(cax)
anim = FuncAnimation(fig,animate,interval = 100, frames = nt)
anim.save('PV.mp4', metadata={'artist':'Guido'},writer='ffmpeg',fps=8,bitrate=500)

plt.show()
plt.close()




