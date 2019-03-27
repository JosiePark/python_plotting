#!/usr/bin/env python

import matplotlib.pyplot as plt

def square_grid_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width):

	sub_width = (1 - left_space - right_space - (ncols-1)*hor_space)/float(ncols)
	sub_height = (1 - top_space - bottom_space - (nrows-1)*ver_space)/float(nrows)

	real_width = fig_width*sub_width
	real_height = real_width # enforce square grid

	height_ratio = sub_height/real_height
	width_ratio = 1./fig_width

	real_top = top_space*height_ratio
	real_bottom = bottom_space*height_ratio

	real_hor = hor_space*fig_width
	real_ver = ver_space*height_ratio

	fig_height = nrows*real_height + real_top + real_bottom + (nrows-1)*real_ver
	
	fig = plt.figure(figsize = (fig_width,fig_height))

	ax = []

	for i in range(ncols):
		for j in range(nrows):
			ax.append(fig.add_axes([left_space+i*sub_width + i*hor_space,bottom_space+(nrows-j-1)*sub_height+(nrows-j-1)*ver_space,sub_width,sub_height]))
		
	return ax

def a4_plot(nrows,ncols,left_space,right_space,bottom_space,top_space,ver_space,hor_space):
# take inputs as coordinates i.e between 0 and 1

	fig_height = 11.69
	fig_width = 8.27

	fig = plt.figure(figsize = (fig_width,fig_height))

	# calculate width and height of plots

	sub_width = (1. - left_space - right_space - (ncols - 1)*hor_space)/float(ncols)
	sub_height = (1. - top_space - bottom_space - (nrows - 1)*ver_space)/float(nrows)

	ax = []

	for i in range(ncols):
		for j in range(nrows):
			ax.append(fig.add_axes([left_space + i*sub_width + i*hor_space,bottom_space + (nrows-j-1)*sub_height + (nrows-j-1)*ver_space,sub_width,sub_height]))
	
	return ax

def spd_grid_plot(left_space,right_space,bottom_space,top_space,ver_space,hor_space):

	fig_height = (11.69/4.)*2.
	fig_width = (8.27/3.)*2.

	nrows = 2
	ncols = 2

	fig = plt.figure(figsize = (fig_width,fig_height))
	
	sub_width = (1. - left_space - right_space - (ncols - 1)*hor_space)/float(ncols)
	sub_height = (1. - top_space - bottom_space - (nrows - 1)*ver_space)/float(nrows)

	ax = []

	for i in range(ncols):
		for j in range(nrows):
			ax.append(fig.add_axes([left_space + i*sub_width + i*hor_space,bottom_space + (nrows-j-1)*sub_height + (nrows-j-1)*ver_space,sub_width,sub_height]))

	return ax
