import numpy as np

def halo_inverse(wt, tab_k, multiscales=False):

	'''
	multiscales allows reconstruction when wavelets have features that
	are not necisacrliy at their original scales
	'''

#----------------definitions----------------#

	ko=5.336
	na=wt.shape[1]
	nb=wt.shape[2]
 
	a=ko/tab_k
	imagetot=np.complex_(np.zeros((na,nb)))


#--------------Coords---------------------------##

	x = np.arange( na , dtype=float )
	y = np.arange( nb , dtype=float )
	x , y = np.meshgrid( x, y )
	if (nb % 2) == 0:
		x = ( 1.*x - (nb)/2. )/ nb 
		shiftx = (nb)/2.
	else:
		x = (1.*x - (nb-1)/2.)/ nb
		shiftx = (nb-1.)/2.+1

	if (na % 2) == 0:
		y = (1.*y-(na/2.))/na
		shifty = (na)/2.
	else:
		y = (1.*y - (na-1)/2.)/ na
		shifty = (na-1.)/2+1

#---------------------transform----------------------#

    
        
	for h in range(tab_k.shape[0]):
		uv = 0
        
		uv = np.exp ( -0.5 * ( abs ( a[h] * np.sqrt ( x**2. + y**2. ) )- ko)**2.)
		uv = uv / a[h]


		imageFT = np.roll((np.fft.fft2(wt[h,:,:])),int(shiftx), axis=1)
		imageFT = np.roll(imageFT, int(shifty), axis=0)
		imageFT = imageFT*uv
		
		imageFT = np.roll(imageFT ,int(shiftx), axis=1)
		imageFT = np.roll(imageFT, int(shifty), axis=0)
		
		image = np.fft.ifft2(imageFT)
		
		imagetot = imagetot+image 

    #imagetot=imagetot*0.95
	return imagetot
