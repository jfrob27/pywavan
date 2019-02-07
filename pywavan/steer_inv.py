import numpy as np
from astropy.io import fits
from .wavan import uv_plane

def steer_inv(steer,Him,Sclim):

	#figfile = '/Users/robitaij/postdoc/Herschel/W43_density_galcut_nan.fits'
	#HDU = fits.open(figfile)
	#im = HDU[0].data
	#Nite = 6
	N=11            #Number of directions (odd)

	#--------------------Definitions----------------------#
	na = im.shape[1]
	nb = im.shape[0]

	x, y, shiftx, shifty, ishiftx, ishifty = uv_plane(na, nb)
	x *= na/2.
	y *= nb/2.
	
	rt = np.sqrt(x**2.+y**2.)
	x[x == 0] = 1
	theta = np.arctan(y/x)

	del x, y

	#-----------------------------------------
	#Angular function
	#-----------------------------------------

	Gk = np.zeros((N,nb,na))
	Gkt = np.zeros((nb,na))

	for j in range(N):

		AA = 2**(N-1) * np.math.factorial(N-1)/np.sqrt(np.float(N*np.math.factorial(2*(N-1))))
	
		#Condition A
		#-------------------
	
		ca = np.where(theta-(np.pi*j)/N < np.pi/2.)
	
		Gkt[ca] = AA * np.cos( theta[ca] - (np.pi*j)/N )**(N-1)

		#Condition B
		#-------------------

		cb = np.where(theta-(np.pi*j)/N >= np.pi/2)

		Gkt[cb] = 0.
	
		#Delete reflection
		#-------------------

		if (j >= 0) & (j <= 2):
			Gkt[0:nb,0:na/2] = 0.

		if (j > 2) and (j < 9):
			Gkt[0:nb/2,0:na] = 0.

		if j >= 9:
			Gkt[0:nb,(na/2)-1:na] = 0.

		#-------------------

		Gk[j,:,:] = Gkt

	Gkt=0

	#-----------------------------------------
	#First division uv-plan
	#-----------------------------------------

	H0 = np.zeros((nb,na))
	L0 = np.zeros((nb,na))

	max = np.pi/2.

	r = rt*max/(na/2.)

	#Condition A
	#-------------------

	ca = np.where((r > np.pi/4.) & (r < np.pi/2.))

	H0[ca] = np.cos( np.pi/2. * np.log(2. * r[ca] / np.pi)/np.log(2.) )

	L0[ca] = 2.*np.cos( np.pi/2. * np.log(4. * r[ca] / np.pi)/np.log(2.) )

	#Condition B
	#-------------------

	cb = np.where(r <= np.pi/4.)

	H0[cb] = 0.

	L0[cb] = 2.

	#Condition C
	#-------------------

	cc = np.where(r >= np.pi/2.)

	H0[cc] = 1.

	L0[cc] = 0.

	#-----------------------------------------
	#Iterative division uv-plan
	#-----------------------------------------

	Hi = np.zeros((Nite,nb,na))
	Li = np.zeros((Nite,nb,na))
	Hit = np.zeros((nb,na))
	Lit = np.zeros((nb,na))

	for ite in range(Nite):

		max = np.pi/2.

		r = (rt*max)/(na/(2.**(ite+1.)))

		#Condition A
		#-------------------

		ca = np.where((r > np.pi/4.) & (r < np.pi/2.))

		Hit[ca] = np.cos( np.pi/2. * np.log(2. * r[ca] / np.pi)/np.log(2.) )

		Lit[ca] = 2.*np.cos( np.pi/2. * np.log(4. * r[ca] / np.pi)/np.log(2.) )

		#Condition B
		#-------------------

		cb = np.where(r <= np.pi/4)

		Hit[cb] = 0.

		Lit[cb] = 2.

		#Condition C
		#-------------------

		cc = np.where(r >= np.pi/2)

		Hit[cc] = 1.

		Lit[cc] = 0.

		#-------------------

		Li[ite,:,:] = Lit
		Hi[ite,:,:] = Hit

		Lit = Lit * 0.
		Hit = Hit * 0.

	#-----------------------------------------
	#Wavelet transform
	#-----------------------------------------

	steer = np.zeros((Nite,N,nb,na), dtype=complex)

	#High-pass filter
	#--------------------

	uvWav = np.copy(H0)

	imFT = np.fft.fft2(im)
	uvplan = np.roll(uvWav,int(ishiftx), axis=1)
	uvplan = np.roll(uvplan,int(ishifty), axis=0)
	WimFT = imFT*np.conj(uvplan)
	Him = np.fft.ifft2(WimFT)

	#Scaling filter
	#--------------------

	Scal = Li[Nite-1,:,:]

	uvWav = Scal/2.

	uvplan = np.roll(uvWav,int(ishiftx), axis=1)
	uvplan = np.roll(uvplan,int(ishifty), axis=0)
	WimFT = imFT * np.conj(uvplan)
	Sclim = np.fft.ifft2(WimFT)

	#Steerable wavelets
	#--------------------

	for ite in range(Nite):

		for j in range(N):

			if ite == 0:
				Bp = Hi[0,:,:] * L0/2.
			else:
				Bp = Hi[ite,:,:] * Li[ite-1,:,:]

			uvWav = Bp * Gk[j,:,:]

			uvplan = np.roll(uvWav,int(ishiftx), axis=1)
			uvplan = np.roll(uvplan,int(ishifty), axis=0)
			WimFT = imFT * np.conj(uvplan)
			steer[ite,j,:,:] = np.fft.ifft2(WimFT)

	return rec	