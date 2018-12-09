import numpy as np
from math import factorial,sqrt,pi

__author__ = "Phi Breen"
__license__ = "MIT"
__email__ = "phil.breen@ed.ac.uk"

def car2sph(x):
	""" Coverts cartesian to spherical coordinates """
	co = np.zeros_like(x)
	co[:,0] = np.sum(x**2,axis=-1)**(1./2)
	co[:,1] = np.arctan2(x[:,1],x[:,0])
	co[:,2] = np.arccos(x[:,2]/co[:,0])
	return co


def alp(l,x,derivative=False):
	""" Calculates all associated Legendre polynomials up to degree l 
	(including all orders up to l, i.e. m<=l) 

	Keyword arguments:
	l -- maximum degree to be calculated
	derivative -- if true also returns diction containing derivative 
	"""
	plm={}
	for mi in range(l+1):
		plmi = 1.0
		if mi > 0:
			plmi = -1.*(2.*mi - 1.)*np.sqrt(1.0 - x**2)*plm[(mi-1,mi-1)]
		plm1m=plmi
		plm2m=0.0
		plm[(mi,mi)] = plmi
		for li in range(mi+1,l+1):
			plmi=(x*(2.*li-1.)*plm1m-(li+mi-1.)*plm2m)/(li-mi)
			plm2m=plm1m
			plm1m=plmi
			plm[(mi,li)] = plmi

	if derivative:
		dplmi = 0.0
		dplm={}
		tc = (x*x-1.0)
		for mi in range(l+1):
			for li in range(mi,l+1):
				if li == 0:
					dplmi = 0.0
				elif mi == li:
					dplmi = li*x*plm[(mi,li)]/tc
				else:
					dplmi = (li*x*plm[(mi,li)]-(li-mi)*plm[(mi,li-1)] )/tc
				dplm[(mi,li)] = dplmi
		return plm, dplm
	return plm



def sph(l,azimuthal,polar,derivative=False):
	""" Calculate spherical harmonics up to degree l including order up to l (i.e. m<=1) 

			
	"""
	sph = {}
	dsphaz = {}
	dsphpol = {}

	if derivative:
		leg, dleg  = alp(l, np.cos(polar), derivative=True)	
	else:
		leg = alp(l, np.cos(polar))
		
	ordterm = []
	for mi in range(l+1):
		ordterm.append(np.exp(azimuthal*mi*1j))

	for li in range(l+1):
		for mi in range(li+1):
			sph[(mi,li)] = sqrt( ((2.*li+1.)/(4.*pi))*(factorial(li-mi)/factorial(li+mi)) )*ordterm[mi]*leg[(mi,li)]
		
			
	if derivative:
		for li in range(l+1):
			for mi in range(li+1):
				dsphaz[(mi,li)] = mi*1j*sph[(mi,li)] 
				dsphpol[(mi,li)] = sqrt( ((2.*li+1.)/(4.*pi))*(factorial(li-mi)/factorial(li+mi)) )*ordterm[mi]*dleg[(mi,li)]*-np.sin(polar)
		return sph,dsphaz,dsphpol	

	return sph




def gegenbauer(n, alpha,x,derivative=False):
	""" Calculates Gegenbauer ultraspherical polynomials  """
	gc=[]
	for ni in range(n+1):			
		if ni == 0:
			gegi = np.ones_like(x)
			gegm2 = gegi
		elif ni == 1:
			gegi = 2.*alpha*x
			gegm1 = gegi				
		else:	
			gegi = (1./ni)*( 2.*(ni+alpha-1.)*x*gegm1  - ( ni+2.*alpha-2.)*gegm2 )
			gegm2 = gegm1
			gegm1 = gegi
		gc.append(gegi)

	if derivative:
		dgc = []
		xf = 1.0/(1-x**2)
		for ni in range(n+1):
			if ni == 0:
				dgegi = np.zeros_like(x)
			elif ni ==1:
				dgegi  = np.ones_like(x)*alpha*2.0
			elif ni ==2:
				dgegi = xf*( -ni*x*gc[ni]  + ( ni+2.*alpha-1.)*gc[ni-1] )
		
			dgc.append(dgegi)
		return gc,dgc
	return gc





