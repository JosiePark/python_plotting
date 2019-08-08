#!/usr/bin/env python
import numpy as np

def autocorrelation(traj,ntau,dt):

	[ndims,npoints,nlayers,nt] = traj.shape
	
	# calculate Lagrangian velocity
	
	u = (traj[:,:,:,1:] - traj[:,:,:,:-1])/dt
	
	# calculate velocity variance
	
	var = np.mean(np.power(u,2),3)
	sigma = np.zeros((ndims,ndims,npoints,nlayers))
	for i in range(ndims):
		for j in range(ndims):
			sigma[i,j,:,:] = np.power(np.multiply(var[i,:,:],var[j,:,:]),2)
			
	# calculate velocity cross covariance
	
	R = np.zeros((ndims,ndims,npoints,nlayers,ntau))
	R_tot = np.zeros((ndims,ndims,nlayers,ntau))
	
	for tau in range(ntau):
		t = nt-tau-1
		tmp = np.zeros((ndims,ndims,npoints,nlayers,t))
		for i in range(ndims):
			for j in range(ndims):
				tmp[i,j,:,:,:] = tmp[i,j,:,:,:] + np.multiply(u[i,:,:,:t],u[j,:,:,tau:t+tau])
			
		R[:,:,:,:,tau] = np.mean(tmp,4)
		R[:,:,:,:,tau] = np.divide(R[:,:,:,:,tau],sigma)
	
	R_tot = np.sum(R,2)/npoints
	

	return R_tot
