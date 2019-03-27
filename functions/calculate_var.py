#!/usr/bin/env python
import numpy as np

def pv_anomaly_from_psi(psi,H1,H2,Rd,basinscale):

	[ii,jj,nlayers] = psi.shape
	
	rel = np.zeros((psi.shape))

	for i in range(ii):
		for j in range(jj):
			im1 = i-1
			ip1 = i+1
			jm1 = j-1
			jp1 = j+1
	
			if im1 == -1:
				im1 = ii-1
			if ip1 == ii:
				ip1 = 0
			if jm1 == -1:
				jm1 = jj-1
			if jp1 == jj:
				jp1 = 0
			for k in range(2):
				rel[i,j,k] = -4*psi[i,j,k] + (psi[im1,j,k] + psi[ip1,j,k] + psi[i,jp1,k] + psi[i,jm1,k])

	uscale = 1.
	scale = basinscale/float(jj)
	
	SS = (scale/Rd)**2.
	S1 = SS/(1+H1/H2)
	S2 = SS - S1

	zeta = np.zeros((psi.shape))

	for i in range(ii):
		for j in range(jj):
			zeta[i,j,0] = rel[i,j,0] - S1*(psi[i,j,0] - psi[i,j,1])
			zeta[i,j,1] = rel[i,j,1] - S2*(psi[i,j,1] - psi[i,j,0])
				

	return zeta

def velocity_from_psi(psi):

	[ii,jj,nlayers] = psi.shape
	u = np.zeros((ii,jj,nlayers))
	v = np.zeros((ii,jj,nlayers))

	for i in range(ii):
		for j in range(jj):
			im1 = i - 1
			ip1 = i + 1
			jm1 = j - 1
			jp1 = j + 1
			if(im1 == -1):
				im1 = ii - 1
			if(ip1 == ii):
				ip1 = 0
			if(jm1 == -1):
				jm1 = jj - 1
			if(jp1 == jj):
				jp1 = 0

			u[i,j,:] = -.5*(psi[i,jp1,:] - psi[i,jm1,:])
			v[i,j,:] = .5*(psi[ip1,j,:] - psi[im1,j,:])
	

	return [u,v]
