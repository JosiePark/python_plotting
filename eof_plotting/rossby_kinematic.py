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
	
home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/EOF/' % regime

eof_file = home_dir + '/STATS/EOF/eof.nc'
eof_data = Dataset(eof_file,'r')

if (regime == 1):
	eof = np.transpose(eof_data.variables['EOFs'][:,:,:,2:4])
else:
	eof = np.transpose(eof_data.variables['EOFs'][:,:,:,1:3])

pc_file = home_dir + '/STATS/EOF/pc.nc'
pc_data = Dataset(pc_file,'r')

if (regime == 1):
	pc = np.transpose(pc_data.variables['PCs'][:,2:4])
else:
	pc = np.transpose(pc_data.variables['PCs'][:,1:3])

# FIND VELOCITY AND PLOT

u = np.zeros((2,ii,jj,2))
v = np.zeros((2,ii,jj,2))

for m in range(2):
	[u[m,:,:,:],v[m,:,:,:]] = velocity_from_psi(eof[m,:,:,:]*pc[m,0])

u_tot = np.zeros((2,ii,jj,2))

u_tot[0,:,:,:] = np.sum(u,axis = 0)
u_tot[1,:,:,] = np.sum(v,axis = 0)

psi = eof[0,:,:,:]*pc[0,0]+eof[1,:,:,:]*pc[1,0]

left_space = .08
top_space = .05
bottom_space = .08
right_space = .05
nrows = 2
ncols = 2
fig_width = 8.27
hor_space = .04
ver_space = .02



xmin,xmax = 0,520
xx,yy = np.mgrid[xmin:xmax:512j,xmin:xmax:512j]



X, Y = np.meshgrid(np.linspace(0, 520., int(ii/8)), np.linspace(0, 520., int(jj/8)))

plt.quiver(X,Y,np.transpose(u_tot[0,0:ii:8,0:jj:8,0]),np.transpose(u_tot[1,0:ii:8,0:jj:8,0]))
fig_name = fig_dir + 'eof_rossby_quiver'
plt.savefig(fig_name)

plt.close()

# FIND AMPLITUDE

start_point = np.zeros((ii))
end_point = np.zeros((jj))

eps = 20.

for i in range(ii):
	j = 0
	while (j < jj and abs(psi[i,j,0]) < 0. + eps):
		j+=1
	start_point[i] = j
	while (j < jj and abs(psi[i,j,0]) > 0. + eps):
		j+=1
	if j == jj-1 or j == jj:
		end_point[i] = 0.
		start_point[i] = 0.
	else:
		end_point[i] = j

amp = max(end_point - start_point)
print(amp)
	

x = np.linspace(0,520.,ii)

plt.scatter(x,start_point,label='start_point')
plt.scatter(x,end_point,label='end_point')
plt.xlim([0,520.])
plt.ylim([0,520.])
plt.legend()
plt.show()
plt.close()

#print('start_point = ',start_point)
#print('end_point = ',end_point)

k = 0

ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

for i in range(ncols):
	for j in range(nrows):
		c = ax[k].pcolor(xx,yy,u_tot[i,:,:,j],cmap = cm.jet)
		if j == 1:
			ax[k].set_xlabel('X (km)')
		else:
			ax[k].set_xticklabels([])
		if i == 0:
			ax[k].set_ylabel('Y (km)')
		else:
			ax[k].set_yticklabels([])

		plt.colorbar(c,ax=ax[k])
		ax[0].set_title('Zonal Velocity')
		ax[2].set_title('Meridional Velocity')

		k+=1
fig_name = fig_dir + 'eof_rossby_velocity'
plt.savefig(fig_name)


plt.close()


# FIND ROTATION SPEED

# FIND PROPAGATION SPEED

# FIND CENTER 

