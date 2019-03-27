#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import square_grid_plot

# INPUT PARAMETERS #

regime = 2

home_dir = '/media/josiepark/Seagate Expansion Drive/PhD/DATA/Saves/%i' % regime
fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/%i/EOF/' % regime

# READ EOF #

pc_file = home_dir + '/STATS/EOF/pc.nc'
pc_data = Dataset(pc_file,'r')

pc = np.transpose(pc_data.variables['PCs'][:])
print(pc.shape)

# PLOT EOF MODES #

left_space = .1
top_space = .02
bottom_space = .04
right_space = .01
nrows = 5
ncols = 2
fig_width = 8.27
hor_space = 0.06
ver_space = 0.02

ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)

#width = (1-left_space-right_space)/float(ncols)
#height = (1-top_space - bottom_space)/float(nrows)

#real_width = fig_width*width
#real_height = real_width # preserve equal figure sizes

#top_ratio = top_space/height
#bottom_ratio = bottom_space/height

time = np.arange(0,20000,1)
#print(time.shape)
	
#fig_height = (2+top_ratio+bottom_ratio)*real_height
#fig = plt.figure(figsize = (fig_width,fig_height))
#ax = []
k = 0
for i in range(ncols):
	for j in range(nrows):
		#ax.append(fig.add_axes([left_space+(i)*width,bottom_space+(1-j)*height,width,height]))
		ax[k].plot(time,pc[:,k],linewidth = 0.2)
		if j == nrows-1:
			ax[k].set_xlabel('Time (days)')
		else:
			ax[k].set_xticklabels([])
		if i == 0:
			ax[k].set_ylabel('PC')
		#else:
		#	ax[k].set_yticklabels([])
		ax[k].ticklabel_format(style = 'sci',axis = 'y',scilimits=(0,0))
		text_x = .5
		text_y = .9
		ax[k].text(text_x,text_y,'Mode %i' % (k+1),horizontalalignment = 'center',verticalalignment = 'center', transform = ax[k].transAxes)

		k+=1
fig_name = fig_dir + 'pc'
plt.savefig(fig_name)
plt.show()

## plot frequency spectrum ## 

time_start = 1000
time_end = 20000

plt.close()
left_space = .1
top_space = .04
bottom_space = .05
right_space = .05
fig_width = 8.27
nrows = 5
ncols = 2
ax = square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,0.02,0.05,fig_width)
k = 0

for i in range(ncols):
	for j in range(nrows):
		#ax.append(fig.add_axes([left_space+(i)*width,bottom_space+(1-j)*height,width,height]))
		ax[k].plot(1./time[time_start:time_end],pc[time_start:time_end,k],linewidth = 0.2)
		if j == nrows-1:
			ax[k].set_xlabel('Frequency (days^(-1))')
		else:
			ax[k].set_xticklabels([])
		if i == 0:
			ax[k].set_ylabel('PC')
		#else:
		#	ax[k].set_yticklabels([])
		ax[k].ticklabel_format(style = 'sci',axis = 'y',scilimits=(0,0))
		text_x = .5
		text_y = .9
		ax[k].text(text_x,text_y,'Mode %i' % (k+1),horizontalalignment = 'center',verticalalignment = 'center', transform = ax[k].transAxes)

		k+=1
fig_name = fig_dir + 'pc_spectrum'
plt.savefig(fig_name)
plt.show()


	
	
	




