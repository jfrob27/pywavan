import numpy as np

def nb_scale(im):
	'''
	Calculate the number of independent scale for a given image
	
	Parameters
	----------
	image: array_like
		Input array, must be 2-dimentional and real
	
	Return
	-------
	Number of scales: scalar
	'''
	ko= 5.336
	delta = (2.*np.sqrt(-2.*np.log(.75)))/ko
	
	nx = np.max(im.shape)
	M=int(np.log(nx)/delta)
	
	return M