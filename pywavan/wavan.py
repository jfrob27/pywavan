import numpy as np
from tqdm.notebook import trange
from .imsmooth import imsmooth
from .edges import apodize, padding, depad
from .uv_plane import uv_plane
from .gauss_segmen import gauss_segmen
import matplotlib.pyplot as plt

def fan_trans(image, reso=1, q=0, qdyn=False, skewl=0.4, zeromean=True, pownorm=True,\
				sigma=None, cutpad=True, smooth=False, angular=False, nan_frame=False, **kwargs):
	'''
	Performs fan transform on 'image' input (Kirby, J. F. (2005),Computers and
	Geosciences, 31(7), 846-864). If an array of spatial scales is not specified
	returns a quasi-orthogonal basis (Robitaille, J.-F. et al. (2014), MNRAS,
	440(3), 2726-2741).
	
	Parameters
	----------
	image : array_like
		Input array, must be 2-dimentional and real
	reso : float, optional
		Resolution of the image in pixel^-1
	q : list, optional
		The dimensionless constant controlling how restrictive is the definition of
		non-Gaussiannities. Typically between 1.8 and 3.0. If set to zero, no
		segmentation is performed by the function. If segmentation is performed, a
		list of q with a number of element equals to the number of scales msut be
		provided.
		Ex.: q=[]
			 q=[2.5]*24
	qdyn : boolean, False by default
		If True, the value of q is minimized according to the skewness of the wavelet
		coefficient distributions as a function of scales.
	skewl : float, 0.4 by default
		Skewness limit value to minimize q.
	zeromean : boolean, True by default
		Subtract image's mean value before performing the wavelet transforms
	pownorm : boolean, True by default
		Normalize the wavelet power spectrum array so thatit can be compared to the
		Fourier power spectrum.
	sigma : variable, optional
		Name of the variable for the uncertainty outputs on the wavelet power spectrum
	cutpad : boolean, True by default
		If True, the function returns the wavelet transformed maps with the dimension
		than the original image. If False, wavelet transformed map sizes correspond to
		the sizes specified in the keyword "arrdim".
	smooth : boolean, False by default
		If True, a smoothing is performed at every scales in order to correct
		segmentation's artefact.
	angular : boolean, False by default
		If True, the returned wavelet coefficients have an angular dependency.
		Ex.: wtcoeff = np.zeros((3*M,N,sy,sx), dtype=complex)
		instead of wtcoeff = np.zeros((3*M,sy,sx), dtype=complex)
	
	Keywords
	--------
	scales : array_like, optional
		Array of spatial scales in terms of Fourier wavenumber k
	apodize : float (0.<radius<1.), optional
		Radius of apodization.
		Ex.: apodize = 0.98, 2% of the edges will be smoothly apodized to zero.
	arrdim : two-element array
		np.array([newy,newx]) sizes including the zero value pixel padding.
	Ndir : integer
		Number of directions to sample (can oversample or undersample for power
		spectrum calculation).
		
	Returns
	-------
	wt : data cube of wavelet coefficients (complex array)
		wt[scales,nx,ny]
		If angular = True -> wt[scales,angles,nx,ny]
		If q != 0 ->
			wtcoeff[0:M,:,:] -> all coefficient
			wtcoeff[M:2*M,:,:] -> non-Gaussian coefficients
			wtcoeff[2*M:3*M,:,:] -> Gaussian coefficients
			where M is the number of scales
	S11a : data cube of <|wt|^2>
	wav_k : Array of spatial scales used for the decomposition
	S1a : Wavelet power spectra
		1-dimensional array -> S11(scales)
		If q !=0 ->
			S1a[0,:] -> total wavelet power spectrum
			S1a[1,:] -> non-Gaussian wavelet power spectrum
			S1a[2,:] -> Gaussian wavelet power spectrum
	q : updated q list if qdyn = True.
	'''
	
	#--------------------Definitions----------------------#
	ko= 5.336
	delta = (2.*np.sqrt(-2.*np.log(.75)))/ko
	na = image.shape[1]
	nb = image.shape[0]
	
	nao = np.copy(na)
	nbo = np.copy(nb)
	
	#-----------------Remove mean value-------------------#
	if zeromean == True:
		image -= np.nanmean(image)
	
	#----------Estimate initial standard deviation------------#
	
	#The initial standard deviation is estimated inorder to calculate the correction
	#factor on the reconstructed map -> C_del = sig_0 / sig_r
	#Constant calculated at the end of this function
	
	sig_0 = np.nanstd(image)
		
	#--------------Apodization--------------------#
	
	if 'apodize' in kwargs:
		tapper = apodize(nb,na,kwargs.get('apodize'))
		image = image *tapper
		
	if 'arrdim' in kwargs:
		arrdim = kwargs.get('arrdim')
		na = arrdim[1]
		nb = arrdim[0]
		image = padding(image,arrdim[0],arrdim[1])
		print(na,nb)
		
	if (nan_frame == True):
		nan_frame_pos = np.where(image != image)
		image = np.nan_to_num(image, copy=False, nan=0)
		
	#--------------Spectral Logarithm--------------------#
	
	if 'scales' in kwargs:
		scales = kwargs.get('scales')
		wav_k = scales
		a2 = 1. / (scales * reso)
		M = scales.size
	
	else:
		nx = np.max(np.array([na,nb]))
		
		M=int(np.log(nx)/delta)
		a2=np.zeros(M)
		a2[0]=np.log(nx)
	
		for i in range(M-1):
			a2[i+1]=a2[i]-delta

		a2=np.exp(a2)
		wav_k = 1. / (a2 * reso)
		
	#-----------------UV-Plane--------------#
	
	x, y, shiftx, shifty, ishiftx, ishifty = uv_plane(na, nb)

	#-----------------Variables--------------#
	if ('Ndir' in kwargs):
		N = kwargs.get('Ndir')  #Number of orientation for the Morlet wavelet (not optimally sampled)
		delta = np.pi/N
	else:
		N = int(np.pi/delta)	#Number of orientation for the Morlet wavelet
	
	if (cutpad == True):
		sx = nao
		sy = nbo
	else:
		sx = na
		sy = nb
		
	if (angular == False):
	
		wt = np.zeros((M,sy,sx), dtype=complex)
		if (q != 0):
			W1c = np.zeros((M,sy,sx), dtype=complex)
			W1n = np.zeros((M,sy,sx), dtype=complex)
			wtcoeff = np.zeros((3*M,sy,sx), dtype=complex)
		else:
			wtcoeff = np.zeros((M,sy,sx), dtype=complex)
			
	else:
	
		wt = np.zeros((M,N,sy,sx), dtype=complex)
		if (q != 0):
			W1c = np.zeros((M,N,sy,sx), dtype=complex)
			W1n = np.zeros((M,N,sy,sx), dtype=complex)
			wtcoeff = np.zeros((3*M,N,sy,sx), dtype=complex)
		else:
			wtcoeff = np.zeros((M,N,sy,sx), dtype=complex)
		
		
	S11 = np.zeros((M,sy,sx))

	if (q != 0):
		S1a = np.zeros((3,M))
		S1c = np.zeros((M,sy,sx))
		S1n = np.zeros((M,sy,sx))
		S11a = np.zeros((3*M,sy,sx))
	else:
		S1a = np.zeros(M)
		S11a = np.zeros((M,sy,sx))
		
	
	
	a = ko * a2				#Scales in the wavelet space
	
	
	#----------------Wavelet transfom------------------------#
	imFT = np.fft.fft2(image)
	#imFTsh = np.fft.fftshift(image)
	imFT= np.roll(imFT,int(shiftx), axis=1)
	imFT= np.roll(imFT,int(shifty), axis=0)
	
	for j in trange(M):
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
			if ('arrdim' in kwargs) & (cutpad == True):
				W1 = depad(W1,nbo,nao)
			
			if (nan_frame == True):
				W1[nan_frame_pos] = float('NaN')
			
			if (angular == False):
				wt[j,:,:]= wt[j,:,:]+ W1
			else:
				wt[j,i,:,:]= wt[j,i,:,:]+ W1
			S11[j,:,:]= S11[j,:,:] + np.abs(W1)**2.
			
	#----------------Segmentation------------------------#
			if q != 0:
				
				cohe, gcoeff, nq = gauss_segmen(W1, q=q[j], qdyn=qdyn, skewl=skewl)
			

				if (W1[gcoeff].shape[0] >  0):
					Wnp = np.zeros((sy,sx), dtype=complex)
					Wnp[gcoeff]=W1[gcoeff]
					if (angular == False):
						W1n[j,:,:] = W1n[j,:,:]+ Wnp
					else:
						W1n[j,i,:,:] = W1n[j,i,:,:]+ Wnp
					S1n[j,:,:] = S1n[j,:,:] + np.abs(Wnp)**2.
					Wnp=Wnp*0.
					
				
				if (W1[cohe].shape[0] > 0):
					Wcp = np.zeros((sy,sx), dtype=complex)
					Wcp[cohe]=W1[cohe]
					if (angular == False):
						W1c[j,:,:] = W1c[j,:,:] + Wcp
					else:
						W1c[j,i,:,:] = W1c[j,i,:,:] + Wcp
					S1c[j,:,:] = S1c[j,:,:] + np.abs(Wcp)**2.
					Wcp=Wcp*0.
		
		if q != 0:
			q[j] = nq
				
		if smooth == True:
			W1ns = W1n[j,:,:]
			W1ns = imsmooth(W1ns, (2.*np.sqrt(2.*np.log(2.)))/(wav_k[j]*reso*2.*np.pi))
			W1n[j,:,:] = W1ns
			del W1ns
			
			W1cs = W1c[j,:,:]
			W1cs = imsmooth(W1cs, (2.*np.sqrt(2.*np.log(2.)))/(wav_k[j]*reso*2.*np.pi))
			W1c[j,:,:] = W1cs
			del W1cs
			
			S1ns = S1n[j,:,:]
			S1ns = imsmooth(S1ns, (2.*np.sqrt(2.*np.log(2.)))/(wav_k[j]*reso*2.*np.pi))
			S1n[j,:,:] = S1ns
			del S1ns
			
			S1cs = S1c[j,:,:]
			S1cs = imsmooth(S1cs, (2.*np.sqrt(2.*np.log(2.)))/(wav_k[j]*reso*2.*np.pi))
			S1c[j,:,:] = S1cs
			del S1cs
				
	#----------------Wavelet power spectra---------------#
		
		if pownorm==True:
			if q != 0:
				#Power spectra with scale power normalisation
				S1a[0,j]=np.nansum(S11[j,:,:]) * a[j]**2. * delta / (float(N) * na * nb)
				S1a[1,j]=np.nansum(S1c[j,:,:]) * a[j]**2. * delta / (float(N) * na * nb)
				S1a[2,j]=np.nansum(S1n[j,:,:]) * a[j]**2. * delta / (float(N) * na * nb)
				
				S11a[j,:,:] = S11[j,:,:] * a[j]**2. * delta / float(N)
				S11a[M+j,:,:] = S1c[j,:,:] * a[j]**2. * delta / float(N)
				S11a[2*M+j,:,:] = S1n[j,:,:] * a[j]**2. * delta / float(N)
				
			else:
				S1a[j]=np.nansum(S11[j,:,:]) * a[j]**2. * delta / (float(N) * na * nb)
				
				S11a[j,:,:] = S11[j,:,:] * a[j]**2. * delta / float(N)
				
		else:
			if q != 0:
				#Power spectra without scale power normalisation
				S1a[0,j]=np.nansum(S11[j,:,:]) * delta / (float(N) * na * nb)
				S1a[1,j]=np.nansum(S1c[j,:,:]) * delta / (float(N) * na * nb)
				S1a[2,j]=np.nansum(S1n[j,:,:]) * delta / (float(N) * na * nb)
				
				S11a[j,:,:] = S11[j,:,:] * delta / float(N)
				S11a[M+j,:,:] = S1c[j,:,:] * delta / float(N)
				S11a[2*M+j,:,:] = S1n[j,:,:] * delta / float(N)
				
			else:
				S1a[j]=np.nansum(S11[j,:,:]) * delta / (float(N) * na * nb)
				
				S11a[j,:,:] = S11[j,:,:] * delta / float(N)
			
	#Calculate uncertainties for the power spectra
	if (sigma == None):
		sigma = None
	else:
		sigma = np.nanstd(S11a, axis=(1,2))
	
	if q != 0:
		if (angular == False):
			#Calculate the correction factor -> C_del = sig_0 / sig_r
			sig_r = np.nanstd(np.sum(wt,axis=0).real)
			C_del = sig_0 / sig_r
			
			wtcoeff[0:M,:,:] = wt * C_del
			wtcoeff[M:2*M,:,:] = W1c * C_del
			wtcoeff[2*M:3*M,:,:] = W1n * C_del
		else:
			sig_r = np.nanstd(np.sum(wt,axis=(0,1)).real)
			C_del = sig_0 / sig_r
			
			wtcoeff[0:M,:,:,:] = wt * C_del
			wtcoeff[M:2*M,:,:,:] = W1c * C_del
			wtcoeff[2*M:3*M,:,:,:] = W1n * C_del
	else:
		if (angular == False):
			sig_r = np.nanstd(np.sum(wt,axis=0).real)
			C_del = sig_0 / sig_r
			wtcoeff = wt * C_del
		else:
			sig_r = np.nanstd(np.sum(wt,axis=(0,1)).real)
			C_del = sig_0 / sig_r
			wtcoeff = wt * C_del
		
	return wtcoeff, S11a, wav_k, S1a, q