import numpy as np

def ind_var(x):
	xb = np.mean(x)
	return x-xb
	
def covar(h,w):
	return np.mean(ind_var(h)*ind_var(w))
