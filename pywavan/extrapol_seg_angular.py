import numpy as np
from .sampling import downsample, rescale
from .noisegen import fbm2d
from .wavan import fan_trans
from .edges import padding, depad
from .gauss_segmen import gauss_segmen
from .imsmooth import imsmooth

def extrapol_seg(coeff, scale, reso, powlawfit, newsize, q=2.5, qdyn=False, skewl=0.4, smooth=False):
	
	fbm = fbm2d(powlawfit[0],newsize, newsize)
	wt, S11a, wav_k, S1a, q = fan_trans(fbm, reso=reso, q=0, qdyn=False, scales=scale)
	del wt, wav_k, q
	
	#Normalisation for the Gaussian extrapolation
	factor = (np.exp(powlawfit[1])*scale**powlawfit[0])/S1a[0]
		
	#Downsampling extrapolation and image wavelet coefficients
	S11abin = downsample(S11a[0,:,:], scale, reso)*factor
	powerbin = downsample(np.abs(coeff)**2., scale, reso)
	
	#Combine image and Gaussian extrapolation and smooth
	powerbinpad = padding(powerbin,S11abin.shape[1],S11abin.shape[0])
	powerbinpad[np.where(powerbinpad == 0.)] = S11abin[np.where(powerbinpad == 0.)]
	powerbinpadsm = rescale(powerbinpad, scale, np.array((newsize,newsize)), reso)
	
	#Perform segmentation
	cohe, gcoeff, nq = gauss_segmen(powerbinpadsm, q=2.5, qdyn=False, skewl=0.4)
	
	coeffpad = padding(coeff,newsize,newsize)

	wtg = np.zeros((newsize,newsize),dtype='complex128')
	wtc = np.zeros((newsize,newsize),dtype='complex128')
	S11ag = np.zeros((newsize,newsize))
	S11ac = np.zeros((newsize,newsize))

	wtc[cohe] = coeffpad[cohe]
	wtg[gcoeff] = coeffpad[gcoeff]
	S11ac[cohe] = powerbinpadsm[cohe]
	S11ag[gcoeff] = powerbinpadsm[gcoeff]

	wtc = depad(wtc,coeff.shape[0],coeff.shape[1])
	wtg = depad(wtg,coeff.shape[0],coeff.shape[1])
	S11ac = depad(S11ac,coeff.shape[0],coeff.shape[1])
	S11ag = depad(S11ag,coeff.shape[0],coeff.shape[1])
	
	if smooth == True:
		wtc = imsmooth(wtc, (2.*np.sqrt(2.*np.log(2.)))/(scale*reso*2.*np.pi))
		wtg = imsmooth(wtg, (2.*np.sqrt(2.*np.log(2.)))/(scale*reso*2.*np.pi))
	
	return wtc, wtg, S11ac, S11ag, powerbinpadsm