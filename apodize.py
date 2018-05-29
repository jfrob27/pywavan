import numpy as np

def apodize(nx, ny, radius):
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

