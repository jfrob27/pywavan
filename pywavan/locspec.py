import numpy as np

def locspec(coords, wav_k, wt, reso, smin, smax, radlim=0.5, phaslim=0.5):
	"""
	Calculate the local wavelet power spectrum for a given position

	Parameters
	----------
	coords : numpy array
		(x, y) pixel coordinate where to calculate the local power spectrum
	wav_k : numpy 1D array
		Array of spatial frequencies (output of fan_trans() function)
	wt : numpy complex array
		data cube of wavelet coefficients with angular dependency
		(output of fan_trans() function with angular option 'True')
	reso : float
		Resolution of the image in pixel^-1
	smin : integer
		lower limit of spatial frequencies (element of the wav_k array)
	smax : integer
		upper limit of spatial frequencies (element of the wav_k array)
		
	Keywords
	--------
	radlim : float, must be < 1.0 and > 0.0
		fraction of radius around 'coords' position for maximum amplitude
		detection --> radius = radlim/(wav_k*reso)
	phaslim: float, must be >= 0.0 and < pi/2
		limit of phase value for maximum amplitude detection
 		

	Returns
	-------

	locpow: numpy 1D array
		the local wavelet power spectrum for wav_k[smin:smax+1]
	xlcoords: numpy 1D array
		x position for maximum amplitude wav_k[smin:smax+1]
	ylcoords: numpy 1D array
		y position for maximum amplitude wav_k[smin:smax+1]
	
	"""
	
	ko= 5.336
	delta = (2.*np.sqrt(-2.*np.log(.75)))/ko
	a = ko/(wav_k*reso)
	Nscl = smax - smin + 1
	Ndir = wt.shape[1]
    
	locpow = np.zeros((Ndir, Nscl))
	xlcoords = np.zeros((Ndir, Nscl))
	ylcoords = np.zeros((Ndir, Nscl))
    
	for i in range(Nscl):
		for j in range(Ndir):
			#Radius around the source (scale/2)
			radius = np.int(radlim/(wav_k[i+smin]*reso))
			#Calculate limits
			#print(coords[1]-radius, coords[0]-radius, coords[1]+radius, coords[0]+radius)
			if (coords[1]-radius <= 0): yl1 = 0
			if (coords[1]-radius > 0): yl1 = coords[1]-radius
			if (coords[0]-radius <= 0): xl1 = 0
			if (coords[0]-radius > 0): xl1 = coords[0]-radius
			if (coords[1]+radius >= wt.shape[2]-1): yl2 = wt.shape[2]
			if (coords[1]+radius < wt.shape[2]-1): yl2 = coords[1]+radius
			if (coords[0]+radius >= wt.shape[3]-1): xl2 = wt.shape[3]
			if (coords[0]+radius < wt.shape[3]-1): xl2 = coords[0]+radius
			#Create sub images for real coefficients and phase
			wtsub = wt[i+smin,j,yl1:yl2, xl1:xl2]
			#wtsub = wt[i+smin,j,coords[1]-radius:coords[1]+radius, coords[0]-radius:coords[0]+radius]
			phasub = np.zeros((yl2-yl1 ,xl2-xl1))
			phasub[wtsub != 0.] = np.abs(np.arctan(wtsub[wtsub != 0.].imag/wtsub[wtsub != 0.].real))
            
			if wtsub.shape == 0:
				print(i,j)
            
			#Locate maximum value in real coefficients inside radius
			maxreal = np.where(wtsub.real == np.max(wtsub.real))
        
			if maxreal[0].size > 1 :
				locpow[j,i] = 'NaN'
			elif (phasub[maxreal] < phaslim) & (maxreal[0] != 0) & (maxreal[1] != 0) \
					& (maxreal[0] != wtsub.shape[0]-1) & (maxreal[1] != wtsub.shape[1]-1):
				locpow[j,i] = np.abs(wtsub[maxreal])**2. * a[i+smin]**2. * delta
				xlcoords[j,i] = maxreal[1]
				ylcoords[j,i] = maxreal[0]
                
	locpow[locpow == 0.] = 'NaN'
    
	return locpow, xlcoords, ylcoords