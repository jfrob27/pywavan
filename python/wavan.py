import numpy as np
from scipy.stats import skew

def gauss_segmen(coeff, q=2.5, qdyn=False, **kwargs):
	
	temoin = np.zeros((coeff.shape[0],coeff.shape[1]))
	if 'index' in kwargs:
		module = np.zeros((coeff.shape[0],coeff.shape[1]))
		module[kwargs.get('index')] = np.abs(coeff[kwargs.get('index')])
	else:
		module = np.abs(coeff)
	tresh = module.max()
	treshp = module.max()*2.

	while ((treshp-tresh) != 0):
		tresh = treshp
		temoin = temoin*0
		
		indx = np.where((module <= tresh) & (module > 0.))
		temoin = (module[indx])**2.
		Sigtresh = np.sum(temoin)/(temoin.shape[0])
		treshp = q *np.sqrt(Sigtresh)
					
		#Adjust q according to the skewness
		if ((treshp-tresh) == 0) & (qdyn==True):
			gcoeff = np.where((module <= tresh) & (module > 0.))
			skewn = skew(np.abs(coeff[gcoeff]))
			if skewn > 0.7:
				q = q - 0.1
				treshp = module.max()*2.
		
	cohe = np.where(module > tresh)
	gcoeff = np.where((module <= tresh) & (module > 0))
			
	return cohe, gcoeff, q

def fan_trans(image, scales=0, reso=1, q=0, qdyn=False, **kwargs):
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
	wav_k : Array of spatial scales used for the decomposition
	S1a : Wavelet power spectrum
		1-dimensional array -> S11(scales)
	'''
	#--------------------Definitions----------------------#
	ko= 5.336
	delta = (2.*np.sqrt(-2.*np.log(.75)))/ko
	na=image.shape[1]
	nb=image.shape[0]
	
	#--------------Spectral Logarithm--------------------#
	
	if scales == 0:
		nx = np.max(image.shape)
		
		M=int(np.log(nx)/delta)
		a2=np.zeros(M)
		a2[0]=np.log(nx)
	
		for i in range(M-1):
			a2[i+1]=a2[i]-delta

		a2=np.exp(a2)
		wav_k = 1. / a2
	else:
		wav_k = scales * reso
		a2 = 1. / scales
		M = scales.size
		
	#-----------------UV-Plane--------------#
	x=np.arange(na)
	y=np.arange(nb)
	x,y=np.meshgrid(x,y)

	if (na % 2) == 0:
		x = (1.*x - (na)/2. )/ na
		shiftx = (na)/2.
		ishiftx = (na)/2.
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
	nS11 = np.zeros((M,nb,na))
	wt = np.zeros((M,nb,na), dtype=complex)

	if (q != 0):
		S1a = np.zeros((3,M))
		S1c = np.zeros((M,nb,na))
		S1n = np.zeros((M,nb,na))
		nS1c = np.zeros((M,nb,na))
		nS1n = np.zeros((M,nb,na))
		W1c = np.zeros((M,nb,na), dtype=complex)
		Wcp = np.zeros((nb,na), dtype=complex)
		W1n = np.zeros((M,nb,na), dtype=complex)
		Wnp = np.zeros((nb,na), dtype=complex)
		S11a = np.zeros((3*M,nb,na))
		wtcoeff = np.zeros((3*M,nb,na), dtype=complex)
	else:
		S1a = np.zeros(M)
		S11a = np.zeros((M,nb,na))
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
			#Wavelet coefficients 
			W1 = np.fft.ifft2(W1FT2)
			
			wt[j,:,:]= wt[j,:,:]+ W1
			nS11[j,:,:]= nS11[j,:,:] + np.abs(W1)
			S11[j,:,:]= S11[j,:,:] + np.abs(W1)**2. 
			
	#----------------Segmentation------------------------#
			if q != 0:
			
				#Set limit based on the noise level of the original
				#power spectrum
		
				cohe, gcoeff, nq = gauss_segmen(W1, q=q[j], qdyn=qdyn)
			

				if (W1[gcoeff].shape[0] >  0):
					Wnp[gcoeff]=W1[gcoeff]
					W1n[j,:,:] = W1n[j,:,:]+ Wnp
					nS1n[j,:,:] = nS1n[j,:,:] + np.abs(Wnp)
					S1n[j,:,:] = S1n[j,:,:] + np.abs(Wnp)**2.
					Wnp=Wnp*0.
					
				if ('double' in kwargs) & (kwargs.get('double') == True) &\
				(W1[cohe].shape[0] > 0):
					#Second iteration of segmentation on non-Gaussiannities\
					#only (coherent part)
					
					cohe2, gcoeff2, coheq = gauss_segmen(W1, q=1.9,\
														 qdyn=False,index=cohe)
					
					Wnp[gcoeff2]=W1[gcoeff2]
					W1n[j,:,:] = W1n[j,:,:]+ Wnp
					nS1n[j,:,:] = nS1n[j,:,:] + np.abs(Wnp)
					S1n[j,:,:] = S1n[j,:,:] + np.abs(Wnp)**2.
					Wnp=Wnp*0.
					
					Wcp[cohe2]=W1[cohe2]
					W1c[j,:,:] = W1c[j,:,:] + Wcp
					nS1c[j,:,:] = nS1c[j,:,:] + np.abs(Wcp)
					S1c[j,:,:] = S1c[j,:,:] + np.abs(Wcp)**2.
					Wcp=Wcp*0.
					
				else:
					if (W1[cohe].shape[0] > 0):
						Wcp[cohe]=W1[cohe]
						W1c[j,:,:] = W1c[j,:,:] + Wcp
						nS1c[j,:,:] = nS1c[j,:,:] + np.abs(Wcp)
						S1c[j,:,:] = S1c[j,:,:] + np.abs(Wcp)**2.
						Wcp=Wcp*0.
					
				q[j] = nq
				
	#----------------Wavelet power spectra---------------#
								
		if q != 0:
            #Power spectra with scale power normalisation
			S1a[0,j]=np.mean(S11[j,:,:]) * a[j]**2. * delta / float(N)
			S1a[1,j]=np.mean(S1c[j,:,:]) * a[j]**2. * delta / float(N)
			S1a[2,j]=np.mean(S1n[j,:,:]) * a[j]**2. * delta / float(N)
            
			S11a[0:M,:,:] = nS11
			S11a[M:2*M,:,:] = nS1c
			S11a[2*M:3*M,:,:] = nS1n
				
		else:
			S1a[j]=np.mean(S11[j,:,:]) * delta / float(N)
			S11a = nS11
			
	if q != 0:
		wtcoeff[0:M,:,:] = wt
		wtcoeff[M:2*M,:,:] = W1c
		wtcoeff[2*M:3*M,:,:] = W1n
	else:
		wtcoeff = wt
		
	return wtcoeff, S11a, wav_k, S1a, q