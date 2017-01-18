import numpy as np

#definition of atrou, fan, halo wavelet functions

###########################################################################


def fan_transform(image, scales=None, reso=None):
	'''
	Performs fan transform on 'image' input (Kirby, J. F. (2005),Computers and
	Geosciences, 31(7), 846–864). If an array of spatial scales is not specified
	returns a quasi-orthogonal basis (Robitaille, J.-F. et al. (2014), MNRAS,
	440(3), 2726–2741). Right now scale option is broken.
	
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
	S11 : Wavelet power spectrum
		1-dimensional array -> S11(scales)
	'''
	#--------------------Definitions----------------------#
	ko= 5.336
    delta = (2.*np.sqrt(-2.*np.log(.75)))/ko
	na=float(image.shape[1])
	nb=float(image.shape[0])
	
	
	#--------------Spectral Logarithm--------------------#
	if reso is None:
		reso = 1
	
	if scales is None:
		nx = np.max(image.shape)
		
		M=int(np.log(nx)/delta)
		a2=np.zeros(M)
		a2[0]=np.log(nx)
	
		for i in range(M-1):
			a2[i+1]=a2[i]-delta

		a2=np.exp(a2)
		tab_k = 1. / a2
	else:
		tab_k	= scales * reso
		
	#-----------------UV-Plane--------------#
	x=np.arange(nb)
	y=np.arange(na)
    x,y=np.meshgrid(x,y)
	if (nb % 2) == 0:
		x=(1.*x - (nb)/2.)/nb
		shiftx = (nb)/2.
	else:
		x= (1.*x - (nb-1)/2.)/nb
		shiftx=(nb-1.)/2.+1

	if (na % 2) == 0:
		y=(1.*y-(na/2.))/na
		shifty=(na)/2.
	else:
		y= (1.*y - (na-1)/2.)/ na
		shifty=(na-1.)/2.+1

	#-----------------Variables--------------#

	S1 = np.zeros((M,nb,na))
	wt = np.zeros((M,nb,na), dtype=complex)
				
    S1a = np.zeros(M)
	
	a = ko*a2				#Scales in the wavelet space
	N = int(np.pi/delta)	#Number of orientation for the Morlet wavelet
	
	
	#----------------Wavelet transfom------------------------#
	imFT = np.fft.fft2(image)
	imFTsh = np.fft.fftshift(image)
	
	for j in range(M):
		for i in range(N):
			uv=0.
			t=float(delta*i)
				
			uv=np.exp( -.5*((a[j]*x - ko*np.cos(t))**2. + (a[j]*y - ko*np.sin(t))**2.))
					
			uv=uv* a[j]		#Energy normalisation
						
			W1FTsh = imFTsh * uv
			W1FT = np.fft.ifftshift(W1FTsh)
			W1 = np.fft.ifft2(W1FT)
			
			wt[j,:,:]= wt[j,:,:]+ W1
			S1[j,:,:]= S1[:,:,j] + np.abs(W1)**2.
										
		S1a[j]=np.mean(S1[j,:,:]) * delta / float(N)

return wt, tab_k, S1a