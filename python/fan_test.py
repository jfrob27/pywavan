import numpy as np

ko= 5.336

na = 256.
nb = 256.

scale = 25.

t = np.pi/4.

#-----------------UV-Plane--------------#
x=np.arange(na)
y=np.arange(nb)
x,y=np.meshgrid(x,y)

if (na % 2) == 0:
	x = (1.*x - (na)/2. )/ na
	shiftx = (na)/2.
else:
	x = (1.*x - (na-1)/2.)/ na
	shiftx = (na-1.)/2.+1
		
if (nb % 2) == 0:
	y = (1.*y-(nb/2.))/nb
	shifty = (nb)/2.
else:
	y = (1.*y - (nb-1)/2.)/ nb
	shifty = (nb-1.)/2+1

#x = np.fft.fftshift(x)
#y = np.fft.fftshift(y)

a = ko * scale

wavFTsh = np.exp( -.5*((a*x - ko*np.cos(t))**2. + (a*y - ko*np.sin(t))**2.))

wavFTsh = wavFTsh * a

wavFT = np.fft.ifftshift(np.complex_(wavFTsh))

#wavFT = np.roll(wavFT, int(shiftx), axis=1)
#wavFT = np.roll(wavFT, int(shifty), axis=0)

wav = np.fft.ifft2(wavFT).real