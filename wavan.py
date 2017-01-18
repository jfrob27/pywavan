import numpy as np

#definition of atrou, fan, halo wavelet functions

###########################################################################


def fan_trans(image, scales=0, reso=1):
	'''
	Performs fan transform on 'image' input (Kirby, J. F. (2005),Computers and
	Geosciences, 31(7), 846-864). If an array of spatial scales is not specified
	returns a quasi-orthogonal basis (Robitaille, J.-F. et al. (2014), MNRAS,
	440(3), 2726-2741). Right now scale option is broken.
	
	Parameters
	----------
	image : array_like
		Input array, must 2-dimentional and real
	scales : array_like, optional
		Array of spatial scales in terms of Fourier wavenumber k
	reso : float, optional
		Resolution of the image in pixel^-1
		
	Returns
	-------
	wt : data cube of wavelet coefficients (complex array)
		wt[scales,nx,ny] -> is the size of the input image
	tab_k : Array of spatial scales used for the decomposition
	S1a : Wavelet power spectrum
		1-dimensional array -> S11(scales)
	'''
	#--------------------Definitions----------------------#
	ko= 5.336
	delta = (2.*np.sqrt(-2.*np.log(.75)))/ko
	na=float(image.shape[1])
	nb=float(image.shape[0])
	
	#--------------Spectral Logarithm--------------------#
	
	if scales == 0:
		nx = np.max(image.shape)
		
		M=int(np.log(nx)/delta)
		a2=np.zeros(M)
		a2[0]=np.log(nx)
	
		for i in range(M-1):
			a2[i+1]=a2[i]-delta

		a2=np.exp(a2)
		tab_k = 1. / a2
	else:
		tab_k = scales * reso
		a2 = 1. / scales
		M = scales.size
		
	#-----------------UV-Plane--------------#
	x=np.arange(na)
	y=np.arange(nb)
	x,y=np.meshgrid(x,y)

	if (na % 2) == 0:
		x = (1.*x - (na)/2. )/ na
		shiftx = (na)/2.
	else:
		x = (1.*x - (na-1)/2.)/ na
		shiftx = (na-1.)/2.+1

	if (nb % 2) == 0:
		y = (1.*y-(nb/2.))/nb
		shifty = (nb)/2.
	else:
		y = (1.*y - (nb-1)/2.)/ nb
		shifty = (nb-1.)/2+1

	#-----------------Variables--------------#

	S11 = np.zeros((M,int(nb),int(na)))
	wt = np.zeros((M,int(nb),int(na)), dtype=complex)

	S1a = np.zeros(M)
	
	a = ko * a2				#Scales in the wavelet space
	N = int(np.pi/delta)	#Number of orientation for the Morlet wavelet
	
	
	#----------------Wavelet transfom------------------------#
	imFT = np.fft.fft2(image)
	#imFTsh = np.fft.fftshift(image)
	imFT= np.roll(imFT,int(shiftx), axis=1)
	imFT= np.roll(imFT,int(shifty), axis=0)
	
	for j in range(M):
		for i in range(N):
			uv=0.
			t=float(delta*i)
				
			uv=np.exp( -.5*((a[j]*x - ko*np.cos(t))**2. + (a[j]*y - ko*np.sin(t))**2.))
					
			uv = uv * a[j]		#Energy normalisation
						
			W1FT = imFT * uv
			W1FT2=np.roll(W1FT,int(shiftx), axis=1)
			W1FT2=np.roll(W1FT2,int(shifty), axis=0)
			W1 = np.fft.ifft2(W1FT2)
			
			wt[j,:,:]= wt[j,:,:]+ W1
			S11[j,:,:]= S11[j,:,:] + np.abs(W1)**2.
										
		S1a[j]=np.mean(S11[j,:,:]) * delta / float(N)

	return wt, tab_k, S1a