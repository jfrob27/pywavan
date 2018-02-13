import numpy as np

def fan_trans(image, scales=0, reso=1, q=0):
	'''
	Performs fan transform on 'image' input (Kirby, J. F. (2005),Computers and
	Geosciences, 31(7), 846-864). If an array of spatial scales is not specified
	returns a quasi-orthogonal basis (Robitaille, J.-F. et al. (2014), MNRAS,
	440(3), 2726-2741). Right now scale option is broken.
	
	Parameters
	----------
	image : array_like
		Input array, must be 2-dimentional and real
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
	na=image.shape[0]
	nb=image.shape[1]
	
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
		ishiftx = (na)/2
	else:
		x = (1.*x - (na-1)/2.)/ na
		shiftx = (na-1.)/2.+1
		ishiftx = (na-1.)/2.

	if (nb % 2) == 0:
		y = (1.*y-(nb/2.))/nb
		shifty = (nb)/2.
		ishifty = (nb)/2
	else:
		y = (1.*y - (nb-1)/2.)/ nb
		shifty = (nb-1.)/2+1
		ishifty = (nb-1.)/2.

	#-----------------Variables--------------#

	S11 = np.zeros((M,nb,na))
	wt = np.zeros((M,nb,na), dtype=complex)

	if q != 0:
		S1a = np.zeros((3,M))
		S1c = np.zeros((M,na,nb))
		S1n = np.zeros((M,na,nb))
		W1c = np.zeros((M,na,nb), dtype=complex)
		Wcp = np.zeros((na,nb), dtype=complex)
		W1n = np.zeros((M,na,nb), dtype=complex)
		Wnp = np.zeros((na,nb), dtype=complex)
		temoin = np.zeros((na,nb))
		module = np.zeros((M,na,nb))
		wtcoeff = np.zeros((3*M,nb,na), dtype=complex)
	else:
		S1a = np.zeros(M)
		wtcoeff = np.zeros((M,nb,na), dtype=complex)
	
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
					
			#uv = uv * a[j]		#Energy normalisation on coefficients directly
						
			W1FT = imFT * uv
			W1FT2=np.roll(W1FT,int(ishiftx), axis=1)
			W1FT2=np.roll(W1FT2,int(ishifty), axis=0)
			W1 = np.fft.ifft2(W1FT2)
			
			wt[j,:,:]= wt[j,:,:]+ W1
			#Wavelet power spectrum with scale power normalisation
			S11[j,:,:]= S11[j,:,:] + np.abs(W1)**2. * a[j]**2. 
			
	#----------------Segmentation------------------------#
	
			if q != 0:
			
				module=abs(W1)
				tresh=module.max()
				treshp=module.max()*2.

				while ((treshp-tresh) != 0):
					tresh=treshp
					temoin = temoin*0
	
					indx=np.where(module <= tresh)
					temoin=(module[indx])**2.
					Sigtresh=np.sum(temoin)/(temoin.shape[0])
					treshp = q *np.sqrt(Sigtresh)
	
				tresh=treshp
				cohe= np.where(module > tresh)

				if (module[cohe].shape[0] > 0):

					Wcp[cohe]=W1[cohe]
					W1c[j,:,:] = W1c[j,:,:] + Wcp
					#Coherent power spectrum with scale power normalisation
					S1c[j,:,:] = S1c[j,:,:] + np.abs(Wcp)**2. * a[j]**2.

					Wcp=Wcp*0
				noncohe =np.where(module <= tresh)

				if (module[noncohe].shape[0] >  0):
					Wnp[noncohe]=W1[noncohe]
					W1n[j,:,:] = W1n[j,:,:]+ Wnp
					#Gaussian power spectrum with scale power normalisation
					S1n[j,:,:] = S1n[j,:,:] + np.abs(Wnp)**2. * a[j]**2.
					Wnp=Wnp*0
				
	#----------------Wavelet power spectra---------------#
								
		if q != 0:
			S1a[0,j]=np.mean(S11[j,:,:]) * delta / float(N)
			S1a[1,j]=np.mean(S1c[j,:,:]) * delta / float(N)
			S1a[2,j]=np.mean(S1n[j,:,:]) * delta / float(N)
		else:
			S1a[j]=np.mean(S11[j,:,:]) * delta / float(N)
			
	if q != 0:
		wtcoeff[0:M,:,:] = wt
		wtcoeff[M:2*M,:,:] = W1c
		wtcoeff[2*M:3*M,:,:] = W1n
	else:
		wtcoeff = wt
		
	return wtcoeff, tab_k, S1a