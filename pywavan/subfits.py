import numpy as np
from astropy.io import fits

def subfits(input,output,coord):
	"""
	Cut a fits image and update the header according to
	the new pixel coordinates

	Parameters
	----------
	input : string
		Path and file name for the fits file to cut
	output : sting
		Path and file name for the new fits file
	coord : numpy array
		4 elements -> np.array([x1,y1,x2,y2])Lower left (x1,y1)
		and upper right (x2,y2) pixels coordnatesfor the new file.
		6 elements -> data cube np.array([x1,y1,z1,x2,y2,z2])
	
	"""
	if coord.shape[0] == 4:
		x1,y1,x2,y2 = coord
	if coord.shape[0] == 6:
		x1,y1,z1,x2,y2,z2 = coord
	if (coord.shape[0] != 4) & (coord.shape[0] != 6):
		print('coord must be a numpy array of 4 or 6 elements')
		return
	
	#Dimension
	width=x2-x1+1
	height=y2-y1+1
	if coord.shape[0] == 6:
		depth=z2-z1+1
		print('Dimensions',width,height,depth)
	else:
		print('Dimensions',width,height)
		
	HDU = fits.open(input)
	im = HDU[0].data
	header = HDU[0].header
	
	header['CRPIX1']-=x1
	header['CRPIX2']-=y1
	if coord.shape[0] == 6:
		header['CRPIX3']-=z1
		
	if coord.shape[0] == 4:
		subim = np.zeros((height,width))
	if coord.shape[0] == 6:
		subim = np.zeros((depth,height,width))
		
	for i in range(width):
		for j in range(height):
			if coord.shape[0] == 4:
				subim[j,i]=im[y1+j,x1+i]
			if coord.shape[0] == 6:
				for k in range(depth):
					subim[k,j,i]=im[z1+k,y1+j,x1+i]
	
	del im
	fits.writeto(output,subim,header,overwrite=True)