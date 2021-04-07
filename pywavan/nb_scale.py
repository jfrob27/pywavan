import numpy as np

def nb_scale(arr):
	'''
	Calculate the number of independent scale for a given image
	
	Parameters
	----------
	dimension: array_like
		Input array, must be 2-elements corresponding to image size
	
	Return
	-------
	Number of scales: scalar
	'''
	ko= 5.336
	delta = (2.*np.sqrt(-2.*np.log(.75)))/ko
	
	nx = np.max(arr)
	M=int(np.log(nx)/delta)
	
	return M