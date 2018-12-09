import numsph
import numpy as np
import scipy.special
from scipy.misc import derivative
import time

def testsph(n,l):
	az = np.random.uniform(0, 2.*np.pi, n)
	pol = np.random.uniform(0, np.pi, n)
	
	t0 = time.time()
	nsph = numsph.sph(l,az,pol)
	t1 = time.time()	
	tot1 = t1-t0	
	
	scisph = {}
	t0 = time.time()
	for li in range(l+1):
		for mi in range(li+1):
			 scisph[(mi,li)] = scipy.special.sph_harm(mi,li,az,pol)
	t1 = time.time()
	tot2 = t1-t0

	err = 0.0
	for ki in nsph.keys():
		err += np.sum(np.abs(nsph[ki]-scisph[ki] ))
	err /= (n*len( nsph.keys()))	

	def sph_harm_az(az):
		return  scipy.special.sph_harm(1,1,az,pol)

	def sph_harm_pol(pol):
		return  scipy.special.sph_harm(1,1,az,pol)
	
	_,dsphaz,dsphpol = numsph.sph(1,az,pol,derivative=True) 
	err2 = np.sum( np.abs(derivative(sph_harm_az,az,dx=1e-6)-dsphaz[(1,1)]))/n
	err3 =  np.sum(np.abs(derivative(sph_harm_pol,pol,dx=1e-6)-dsphpol[(1,1)]))/n
	return tot1,tot2,err,err2,err3


def testgeg(n,alpha,nx):
	x = np.random.uniform(-1.0, 1.0, nx)

	t0 = time.time()
	geg1 = numsph.gegenbauer(n,alpha,x)
	t1 = time.time()	
	tot1 = t1-t0	
	
	scigeg = []
	t0 = time.time()
	for ni in range(n+1):
		gegi = scipy.special.gegenbauer(ni,alpha)
		scigeg.append(gegi(x))
	t1 = time.time()
	tot2 = t1-t0

	err = 0.0
	for ki in range(len(gegi)):
		err += np.sum(np.abs(geg1[ki]-scigeg[ki] ))
	err /= (nx*len(gegi))	
	geg1,dgeg1 =  numsph.gegenbauer(n,alpha,x,derivative=True)
	err2 = np.sum(np.abs(derivative(gegi,x,dx=0.01)-dgeg1[-1]))/nx
	return tot1,tot2,err,err2


def testall():
	print("# testing sph ")
	print("# numsphtime scipytime diff diffaz diffpol n speedup")
	for n in range(3,7):
		t1, t2,err,err2,err3 = testsph(10**n,6)
		print(" {:3.2e} {:3.2e} {:3.2e} {:3.2e} {:3.2e} {} {:3.2f} ".format(t1,t2,err,err2,err3,10**n,t2/t1))
	print("# gegenbauer ")
	print("# numsphtime scipytime diff diff_dev  n speedup")	
	for n in range(0,6):
		t1, t2,err,derr = testgeg(1,2.0,10**n)
		print(" {:3.2e} {:3.2e} {:3.2e} {:3.2e} {} {:3.2f} ".format(t1,t2,err,derr,10**n,t2/t1))

if __name__ == '__main__':
	testall()



