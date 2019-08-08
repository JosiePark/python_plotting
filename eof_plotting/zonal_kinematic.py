#!/usr/bin/env python

## CODE TO HELP CONSTRUCT ZONAL EOF ANALYTICALLY

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

regime = 2
ii = 512
jj = 512

# Read the zonal eof pattern from file

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/EOF/' % regime

eof_file = home_dir + '/STATS/EOF/eof.nc'
eof_data = Dataset(eof_file,'r')
eof = np.transpose(eof_data.variables['EOFs'][:,:,:,0])

xmin,xmax = 0,520
xx,yy = np.mgrid[xmin:xmax:512j,xmin:xmax:512j]

plt.pcolor(xx,yy,eof[:,:,0],cmap=cm.jet)
plt.colorbar()
plt.show()
plt.close()

plt.plot(eof[0,:,0],np.linspace(0,ii,ii))
plt.show()
plt.close()

## find coordinate of minimum and maximum eof value

idMin = np.argmin(eof[0,:,0])
idMax = np.argmax(eof[0,:,0])

y = np.linspace(0,ii,ii)

print(idMin,idMax)

p3 = np.poly1d(np.polyfit(y,eof[0,:,0],3))
p4 = np.poly1d(np.polyfit(y,eof[0,:,0],4))

plt.plot(y,eof[0,:,0],'.',y,p3(y),'--',y,p4(y),'-')
plt.legend(['EOF','p3','p4'])
plt.show()
plt.close()

print(np.polyfit(y,eof[0,:,0],3))
