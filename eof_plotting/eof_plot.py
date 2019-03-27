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
ii = 512
jj = ii

# READ EOF #

eof_file = home_dir + '/STATS/EOF/eof.nc'
eof_data = Dataset(eof_file,'r')
eof = np.transpose(eof_data.variables['EOFs'][:])
print(eof.shape)

# PLOT EOF MODES #

left_space = .08
top_space = .05
bottom_space = .15
right_space = .05
nrows = 2
ncols = 5
fig_width = 8.27

#ax = square_grid_plot(nrows,ncold,left_space,right_space,bottom_space,top_space,0.,0.,fig_width)

width = (1-left_space-right_space)/float(ncols)
height = (1-top_space - bottom_space)/float(nrows)

real_width = fig_width*width
real_height = real_width # preserve equal figure sizes

top_ratio = top_space/height
bottom_ratio = bottom_space/height

fig_height = (2+top_ratio+bottom_ratio)*real_height


xmin,xmax = 0,520
xx,yy = np.mgrid[xmin:xmax:512j,xmin:xmax:512j]

for l in range(2):
	fig = plt.figure(figsize = (fig_width,fig_height))
	ax = []
	k = 0
	for i in range(5):
		for j in range(2):

			ax.append(fig.add_axes([left_space+(i)*width,bottom_space+(1-j)*height,width,height]))
			ax[k].pcolor(xx,yy,eof[k,:,:,l],cmap = cm.jet)
			if j == 1:
				ax[k].set_xlabel('X (km)')
			else:
				ax[k].set_xticklabels([])
			if i == 0:
				ax[k].set_ylabel('Y (km)')
			else:
				ax[k].set_yticklabels([])
			text_x = .5
			text_y = .9
			ax[k].text(text_x,text_y,'Mode %i' % (k+1),horizontalalignment = 'center',verticalalignment = 'center', transform = ax[k].transAxes)
			k+=1
	
	if (l==0):
		fig_name = fig_dir + 'top_eof'
	else:
		fig_name = fig_dir + 'bottom_eof'
	plt.savefig(fig_name)
	plt.show()
	plt.close(fig)


	
	
	




