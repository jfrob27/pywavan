import numpy as np
from scipy.stats import skew

def gauss_segmen(coeff, q=2.5, qdyn=False, skewl=0.4):
	
	module = np.abs(coeff)
	tresh = module.max()
	treshp = module.max()*2.

	while ((treshp-tresh) != 0):
		tresh = treshp
		
		indx = np.where((module <= tresh) & (module > 0.))
		Sigtresh = np.mean((module[indx])**2.)
		treshp = q *np.sqrt(Sigtresh)
					
		#Adjust q according to the skewness
		if ((treshp-tresh) == 0) & (qdyn==True):
			gcoeff = np.where((module <= tresh) & (module > 0.))
			skewn = skew(np.abs(coeff[gcoeff]))
			if skewn > skewl:
				q = q - 0.1
				treshp = module.max()*2.
		
	cohe = np.where(module > tresh)
	gcoeff = np.where((module <= tresh) & (module > 0))
			
	return cohe, gcoeff, q