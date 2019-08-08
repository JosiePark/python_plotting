#!/usr/bin/env python

## PLOTS STOKES DRIFT PROFILE ##

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import sys
sys.path.insert(0,'/home/josiepark/Project/PhD/CODE/PYTHON_PLOTTING/functions/')
from grid_plot import square_grid_plot

ii = 512
jj = 512

fig_dir = '/home/josiepark/Project/PhD/PYTHON_FIGURES/COMBO/KINEMATIC/'

factor = [[3.*np.pi,6.54],[7.39394,5.4747]]
A = [[526.1,146.27],[609.8531,191.97748]]
T = [48.,44.]
y_c = [302.,282.]
	
k = 2.
lambd = float(ii)/k
c = [lambd/i for i in T]
cnd = [C*2*np.pi/float(ii) for C in c]

func = np.zeros((2,2,512))

y = np.linspace(0,511,512)
for i in range(2):
	for j in range(2):
		f = lambda x : A[i][j]*(1./np.cosh(factor[i][j]*(x-y_c[i])/float(ii)))**2
		b = factor[i][j]/float(ii)
		fyy = lambda x : -(A[i][j]*b*2/np.cosh(b*(x-y_c[i]))**2)**2*(1./(np.cosh(b*(x-y_c[i]))**2) - 4*np.tanh(b*(x-y_c[i]))**2)/(4*cnd[i])
		func[i,j,:] = fyy(y)

left_space = 0.08
right_space = 0.05
bottom_space = 0.08
top_space = 0.05
hor_space = 0.05
ver_space = 0.05
fig_width = 8.27

ax = square_grid_plot(2,2,left_space,right_space,bottom_space,top_space,hor_space,ver_space,fig_width)
k = 0

for i in range(2):
	for j in range(2):
		ax[k].plot(func[i,j,:],y,'b-')
		ax[k].grid()
		ax[k].set_xlim([-1500,1000])
		k+=1	
ax[0].set_title('Coherent Jet')
ax[2].set_title('Latent Jet')
ax[0].set_ylabel('Top Layer')
ax[1].set_ylabel('Bottom Layer')
ax[1].set_xlabel('Stokes Drift')
ax[3].set_xlabel('Stokes Drift')
			
fig_name = fig_dir + 'stokes_drift.png'
plt.savefig(fig_name)
plt.show()
plt.close()



