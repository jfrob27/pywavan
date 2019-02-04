import numpy as np

def apodize(ny, nx, radius):
	"""
	Create edges apodization tapper

	Parameters
	----------
	nx, ny : integers
		size of the tapper
	radius : float
		radius must be lower than 1 and greater than 0.

	Returns
	-------

	tapper : numpy array ready to multiply on your image
				to apodize edges
	"""
	
	if (radius >= 1) or (radius <= 0.):
		print('Error: radius must be lower than 1 and greater than 0.')
		return
	
	ni = np.fix(radius*nx)
	dni = int(nx-ni)
	nj = np.fix(radius*ny)
	dnj = int(ny-nj)
	
	tap1d_x = np.ones(nx)
	tap1d_y = np.ones(ny)
	
	tap1d_x[0:dni] = (np.cos(3. * np.pi/2. + np.pi/2.* (1.* np.arange(dni)/(dni-1)) ))
	tap1d_x[nx-dni:] = (np.cos(0. + np.pi/2. * (1.* np.arange(dni)/(dni-1)) ))
	tap1d_y[0:dnj] = (np.cos(3. * np.pi/2. + np.pi/2. * (1.* np.arange( dnj )/(dnj-1)) ))
	tap1d_y[ny-dnj:] = (np.cos(0. + np.pi/2. * (1.* np.arange(dnj)/(dnj-1)) ))
	
	tapper = np.zeros((ny, nx))
	
	for i in range(nx):
		tapper[:,i] = tap1d_y
		
	for i in range(ny):
		tapper[i,:] = tapper[i,:] * tap1d_x

	return tapper
	
###############################################
	
def padding(input, y, x):

	width = input.shape[1]
	height = input.shape[0]

	output = np.zeros((y,x))
	
	xpos = np.int(x/2 - width/2)
	ypos = np.int(y/2 - height/2)
	
	output[ypos:height+ypos,xpos:width+xpos] = input
	
	return output
	
###############################################

def depad(input, y, x):

	width = input.shape[1]
	height = input.shape[0]
	
	output = np.zeros((y,x))
	
	xpos = np.int(width/2 - x/2)
	ypos = np.int(height/2 - y/2)
	
	output = input[ypos:y+ypos,xpos:x+xpos]
	
	return output