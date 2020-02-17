import numpy as np
from .rebin import congrid
from .uv_plane import uv_plane

def downsample(coeff, scl, reso):
	
	ko= 5.336
	nxbin = np.round(ko * coeff.shape[1] * scl *reso)
	nybin = np.round(ko * coeff.shape[0] * scl *reso)
	
	coefbin = congrid(coeff,(nybin,nxbin),method='linear',centre=True)
	
	return coefbin

#################################################

def rescale(coefbin, scl, newsize, reso):
	
	ko= 5.336
	ratio0 = newsize[0]/coefbin.shape[0]+1
	ratio1 = newsize[1]/coefbin.shape[1]+1
	
	#Resize the pixelated image with more pixels with function repeat
	#Assign a specific number of pixels with congrid
	coefrsz = congrid(np.repeat(np.repeat(coefbin,ratio0,axis=0),ratio1,axis=1),newsize,method='linear',centre=True)
	#coefrsz = np.repeat(np.repeat(coefbin,ratio0,axis=0),ratio1,axis=1)
	
	# Shift to counter part
	coefrsz = np.roll(coefrsz,-np.int(ratio1/4), axis=1)
	coefrsz = np.roll(coefrsz,-np.int(ratio0/4), axis=0)
	
	#Smooth image with the Gaussian envelop of the wavelet
	a = 1. / (scl*reso)
	u, v, shiftx, shifty, ishiftx, ishifty = uv_plane(coefrsz.shape[1], coefrsz.shape[0])
	r = np.sqrt(u**2. + v**2.)
	
	kernel = (np.pi/2.)*np.exp(-.5*(r*a)**2.)
	
	
	coefrszFT = np.fft.fft2(coefrsz)
	coefrszFT= np.roll(coefrszFT,int(shiftx), axis=1)
	coefrszFT= np.roll(coefrszFT,int(shifty), axis=0)
	
	smoothFT = coefrszFT * kernel
	
	smoothFT2=np.roll(smoothFT,int(ishiftx), axis=1)
	smoothFT2=np.roll(smoothFT2,int(ishifty), axis=0) 
	smooth = np.fft.ifft2(smoothFT2)
	
	return smooth.real