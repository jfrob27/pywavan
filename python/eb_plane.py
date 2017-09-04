import numpy as np

def eb_plane(Q,U):

	na=float(Q.shape[1])
	nb=float(Q.shape[0])

	x=np.arange(na)
	y=np.arange(nb)

	x,y=np.meshgrid(x,y)

	if (nb % 2) == 0:
		x = (1.*x - (na)/2.)
	else:
		x = (1.*x - (na-1)/2.)
			
	if (na % 2) == 0:
		y = (1.*y-(nb/2.))
	else:
		y= (1.*y - (nb-1)/2.)

	phi = (-1.)*np.arctan2(y,x) 	# (-1.) factor is coming from the comment following eq. (1) in the
									# Seljak 1997 paper. \phi_{k,n} can be approximated with -1.*\phi_{k}
									# around the pole (or in the flat sky limit). However, I am not sure why it
									# does not appear again in eq. (8).

	#----------------Fourier Transform-------------------#

	QFT = np.fft.fft2(Q)
	QFTsh = np.fft.fftshift(QFT)

	UFT = np.fft.fft2(U)		# (-1.) factor to respect the IAU convention.
	UFTsh = np.fft.fftshift(UFT)

	#-----------------E/B composition------------------#

	#See equations 24 & 25 of Seljak 1997, ApJ, 482, 6S

	angle = 2.*phi #+ np.pi/2.		#For some obscure reasons, + np.pi/2. is needed

	EFTsh = QFTsh*np.cos(angle) + UFTsh*np.sin(angle)

	BFTsh = -1.*QFTsh*np.sin(angle) + UFTsh*np.cos(angle)


	EFT = np.fft.ifftshift(EFTsh)
	BFT = np.fft.ifftshift(BFTsh)

	E = np.fft.ifft2(EFT).real*(-1.)
	B = np.fft.ifft2(BFT).real*(-1.)

	return E, B


