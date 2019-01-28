import numpy as np

def gauss_beam(sigma, ny, nx, FWHM=False):
	'''
	Return a Gaussian beam at the center of an image
	
	Parameters
	----------
	sigma: float 
		Standard deviation (sigma) of the Gaussian function
	nx, ny: integers
		Size of the output image
	FWHM: bool
		If True, the sigma parameter becomes the Full Width Half Max
		
	Return
	----------
	Gaussian beam at the center of the image
	'''
	
	X=np.arange(nx)
	Y=np.arange(ny)
	ymap,xmap=np.meshgrid(X,Y)
	
	if (nx % 2) == 0:
		xmap = xmap - (nx)/2.
	else:
		xmap = xmap - (nx-1.)/2.

	if (ny % 2) == 0:
		ymap = ymap - (ny)/2.
	else:
		ymap = ymap - (ny-1.)/2.
	
	map = np.sqrt(xmap**2.+ymap**2.)
	
	if FWHM == True:
	
		sigma = sigma / (2.*np.sqrt(2.*np.log(2.)))
	
	gauss = np.exp(-0.5*(map)**2./sigma**2.)
		
	return gauss