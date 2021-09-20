import numpy as np
from scipy.stats import skew

def gauss_segmen(coeff, q=2.5, qdyn=False, skewl=0.4):
	'''
	Non- Gaussian segmentation algorythm based on
	Azzalini, A., Farge, M., & Schneider, K. 2005
	Applied and Computational Harmonic Analysis, 18 (Elsevier), 177
	'''
	
	module = np.abs(coeff)
	tresh = np.nanmax(module)
	treshp = np.nanmax(module)*5.

	while (((treshp-tresh) != 0) & (q > 1.)):
		tresh = treshp
		
		indx = np.where((module <= tresh) & (module > 0.))
		Sigtresh = np.nanmean((module[indx])**2.)
		treshp = q *np.sqrt(Sigtresh)
					
		#Adjust q according to the skewness
		if ((treshp-tresh) == 0) & (qdyn==True):
			gcoeff = np.where((module <= tresh) & (module > 0.))
			skewn = skew(np.abs(coeff[gcoeff]))
			if skewn > skewl:
				q = q - 0.1
				treshp = np.nanmax(module)*2.
		
	cohe = np.where(module > tresh)
	gcoeff = np.where((module <= tresh) & (module > 0))
			
	return cohe, gcoeff, q