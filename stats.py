import numpy as np

def ind_var(x);
	xb = np.mean(x)
	return [y- xb for y in x]
	
def covar(h,w):
	n = len(h)
	return np.dot(ind_var(h),ind_var(w)/(n-1)