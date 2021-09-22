import numpy as np
from astropy.io import fits

def nan_frame_cut(input,output):
	'''
	Remove the "NaN frame" of an image and replace edge irregularities
	with zeros. Corrections on the headers are done.
	
	Parameters
	----------
	input : string
		Path and file name for the fits file to cut
	output : sting
		Path and file name for the new fits file
	'''
	HDU = fits.open(input)
	image = HDU[0].data
	header = HDU[0].header
	
	none_nans = np.where(image == image)
	x1,y1 = np.min(none_nans[1]),np.min(none_nans[0])
	x2,y2 = np.max(none_nans[1]),np.max(none_nans[0])
	
	header['CRPIX1']-=x1
	header['CRPIX2']-=y1

	imcrop = np.zeros((y2-y1+1,x2-x1+1))
	imcrop = image[y1:y2+1,x1:x2+1]
	
	del image
	
	fits.writeto(output,imcrop,header,overwrite=True)