import numpy as np
from scipy.optimize import leastsq

def gaussian(x,c,mu,sig):
	'''
	Return Gaussian function
	
	Parameters
	----------
	x : array_like
		x-axis
	c : float
		Scale of the Gaussian
	mu : float
		Mean value of the Gaussian
	sig : float
		Standard deviation of the Gaussian
	'''
	func = c * np.exp( - (x - mu)**2.0 / (2.0 * sig**2.0) )
	return func

def multi_gaussians(x, params):
	'''
	Return the sum of multiple Gaussians according to the parameters
	
	Parameters
	----------
	x : array_like
		x-axis
	params : list
		List of Gaussian parameters : scale, mean, sigma
		Ex.: 	params = []
				params.append([10.0, -0.82, 10.])
				params.append([1.0, -11.54, 1.65])
				params.append([32.0, -2.47, 1.65])
				params.append([60.0, 1.65, 2.47])
				params.append([45.0, 3.0, 2.5])	
	'''
	nb = len(params) / 3
	res = np.zeros(np.size(x))
	for j in range(nb):
		c = params[3*j]
		mu = params[3*j+1]
		sig = params[3*j+2]
		res += gaussian(x,c,mu,sig)
	return res
		
def gaussian_fit( params ):
    fit = multi_gaussians( x, params )
    return (fit - y_proc) 
    
def gaussfit(tab, spec, params):
	'''
	Perform multiple Gaussian fit on a spectrum
	
	Parameters
	----------
	tab : array_like
		Spectrum x-axis, 1D array
	spec : array_like
		Spectrum data, 1D array
	params : list
		List of Gaussian parameters : scale, mean, sigma
		Ex.: 	params = []
				params.append([10.0, -0.82, 10.])
				params.append([1.0, -11.54, 1.65])
				params.append([32.0, -2.47, 1.65])
				params.append([60.0, 1.65, 2.47])
				params.append([45.0, 3.0, 2.5])
	
	Returns
	-------
	fit : tuple -> list, scalar
		Optimized parameters for the Gaussians
		
	History
	-------
	J.-F. Robitaille September 2017
	'''
	global x
	global y
	global y_proc
	x = np.copy(tab)
	y = np.copy(spec)
	y_proc = np.copy(y)
	y_proc[y_proc < 0.1] = 0.0
	fit = leastsq( gaussian_fit, params)
	return fit