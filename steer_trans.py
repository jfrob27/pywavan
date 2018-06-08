import numpy as np
from astropy.io import fits

#def steer_trans(im,Nite=2):

figfile = '/Users/robitaij/postdoc/Herschel/W43_density_galcut_nan.fits'
HDU = fits.open(figfile)
im = HDU[0].data
Nite = 6

#--------------------Definitions----------------------#
na = im.shape[1]
nb = im.shape[0]

x=np.arange(na)
y=np.arange(nb)
x,y=np.meshgrid(x,y)

#-----------------UV-Plane--------------#

if (na % 2) == 0:
	x = (1.*x - (na)/2. )
	shiftx = (na)/2.
	ishiftx = (na)/2.
else:
	x = (1.*x - (na-1)/2.)
	shiftx = (na-1.)/2.+1
	ishiftx = (na-1.)/2.

if (nb % 2) == 0:
	y = (1.*y-(nb/2.))
	shifty = (nb)/2.
	ishifty = (nb)/2
else:
	y = (1.*y - (nb-1)/2.)
	shifty = (nb-1.)/2+1
	ishifty = (nb-1.)/2.
	
rt = np.sqrt(x**2.+y**2.)
x[x == 0] = 1
theta = np.arctan(y/x)

del x, y

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

Hi = np.zeros((nb,na,Nite))
Li = np.zeros((nb,na,Nite))
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

	Li[:,:,ite] = Lit
	Hi[:,:,ite] = Hit

	Lit = Lit * 0.
	Hit = Hit * 0.

#-----------------------------------------
#Wavelet transform
#-----------------------------------------

steer = np.zeros((nb,na,Nite), dtype=complex)

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

Scal = Li[:,:,Nite-1]

uvWav = Scal/2.

uvplan = np.roll(uvWav,int(ishiftx), axis=1)
uvplan = np.roll(uvplan,int(ishifty), axis=0)
WimFT = imFT*np.conj(uvplan)
Sclim = np.fft.ifft2(WimFT)

#Steerable wavelets
#--------------------

for ite in range(Nite):

	#for j=0,N-1 do begin

	if ite == 0:
		Bp = Hi[:,:,0] * L0/2.
	else:
		Bp = Hi[:,:,ite] * Li[:,:,ite-1]

	uvWav=Bp #*Gk[*,*,j]

	uvplan = np.roll(uvWav,int(ishiftx), axis=1)
	uvplan = np.roll(uvplan,int(ishifty), axis=0)
	WimFT = imFT*np.conj(uvplan)
	steer[:,:,ite] = np.fft.ifft2(WimFT)

	#return steer,Him,Sclim	