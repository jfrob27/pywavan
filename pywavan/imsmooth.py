import numpy as np

def imsmooth(im, FWHM):
	'''
	Return the smoothed image
	
	Parameters
	----------
	im: numpy array (2D)
		Image to smooth
	FWHM: float
		Full Width Half Max for the Gaussian function
	'''

	na=float(im.shape[1])
	nb=float(im.shape[0])

	sigma = FWHM / (2.*np.sqrt(2.*np.log(2.)))

	#-----------------UV-Plane--------------#
	X=np.arange(na)
	Y=np.arange(nb)
	u,v=np.meshgrid(X,Y)

	if (na % 2) == 0:
		u = (u - (na)/2. )/na
	else:
		u = (u - (na-1)/2.)/na

	if (nb % 2) == 0:
		v = (v - (nb/2.))/nb
	else:
		v = (v - (nb-1)/2.)/nb
	
	uvmap = np.sqrt(u**2.+v**2.)

	FTgauss = np.pi/2. * np.exp(-(np.pi*uvmap)**2.*2.*sigma**2.)

	#-----------------Convolution-----------#

	FTsmooth = np.fft.fft2(im) * np.fft.ifftshift(FTgauss)

	smooth = np.fft.ifft2(FTsmooth)

	return np.real(smooth)