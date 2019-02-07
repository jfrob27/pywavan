import numpy as np
from .wavan import uv_plane

def halo_transform(image):
	'''
	Performs halo wavelet transform on image.
	Returns wavelets wt as image cube
	'''

	na = image.shape[1]
	nb = image.shape[0]
	ko = 5.336
	delta = (2.*np.sqrt(-2.*np.log(.75)))/ko

	x, y, shiftx, shifty, ishiftx, ishifty = uv_plane(na, nb)

	#--------------Spectral Logarithm--------------------#
	nx = np.max(image.shape)
	M = int(np.log(nx)/delta)
	a2 = np.zeros(M)
	a2[0] = np.log(nx)

	for i in range(M-1):
		a2[i+1] = a2[i]-delta

	a2 = np.exp(a2)
	tab_k = 1. / (a2)
	wt = np.zeros((M,nb,na))

	a = ko*a2

	imageFT = np.fft.fft2(image)
	imageFT= np.roll(imageFT,int(shiftx), axis=1)
	imageFT= np.roll(imageFT,int(shifty), axis=0)

	for j in range(M):
		uv = 0

		uv = np.exp( -0.5*( a[j]*np.sqrt(x**2.+y**2.) - ko )**2. )

		#uv = uv * a[j]

		W1FT = imageFT*(uv)
		W1F2 = np.roll(W1FT,int(ishiftx), axis =1)
		W1F2 = np.roll(W1F2,int(ishifty),axis=0)
		W1 = np.fft.ifft2(W1F2)
		wt[j,:,:] += W1.real

	return wt, tab_k