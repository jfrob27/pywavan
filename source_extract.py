from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from wavan import uv_plane
from halo_transform import halo_transform

figfile = '/Users/robitaij/postdoc/Herschel/Ken_Marsh/CMa_R1_cut_nan.fits'

#Read map
HDU = fits.open(figfile)
im = HDU[0].data
header = HDU[0].header

model = np.zeros(im.shape)
modmap = np.zeros(im.shape)

x, y, shiftx, shifty, ishiftx, ishifty = uv_plane(im.shape[1], im.shape[0])

#Wavelet transform of the image at a specific scale
ko = 5.336
wt,tab_k =halo_transform(im-np.mean(im))
a = ko/tab_k
M = tab_k.shape[0]
wtalt = np.zeros((M,im.shape[0],im.shape[1]))

start = 4
end = 24
wtalt[0:start,:,:] = wt[0:start,:,:]
wtalt[end:M,:,:] = wt[end:M,:,:]

for i in range(start,end):	

	Walt = np.copy(wt[i,:,:])

	valmax = np.max(Walt)
	indmax = np.where(Walt == valmax)

	while Walt[indmax]/np.std(wt[i,:,:]) > 5.:
		
		model[indmax] = 1.
		
		wav = np.exp( -0.5*( a[i]*np.sqrt(x**2.+y**2.) - ko )**2. )
		
		modelFT = np.fft.fft2(model)
		modelFT= np.roll(modelFT,int(shiftx), axis=1)
		modelFT= np.roll(modelFT,int(shifty), axis=0)
		WmodFT = modelFT * wav
		WmodFT2=np.roll(WmodFT,int(ishiftx), axis=1)
		WmodFT2=np.roll(WmodFT2,int(ishifty), axis=0)
		#Wavelet coefficients 
		Wmod = np.fft.ifft2(WmodFT2).real
		modmap += Wmod
		#Model normalisation (not sure why I have to do this)
		Wmod = Wmod/Wmod[indmax] * Walt[indmax]
	
		Walt -= Wmod
		model *= 0.
		valmax = np.max(Walt)
		indmax = np.where(Walt == valmax)
		
	wtalt[i,:,:] = Walt
	
#Image reconstruction
rec = np.sum(wtalt,axis=0)

#Plot

#histo, edges = np.histogram(Walt[Walt != 0.],bins=250)

#plt.figure(figsize=(7,7))
#plt.bar(edges[:-1], histo, width=np.diff(edges), align="edge",edgecolor='none',facecolor='#858bc7',label='PDF')

#plt.figure(figsize=(15,8))

#plt.subplot(121)
#plt.imshow(Wim,origin='below',cmap='gray')
#plt.subplot(122)
#plt.imshow(Walt,origin='below',cmap='gray')

fits.writeto("/Users/robitaij/postdoc/Herschel/Ken_Marsh/analysis/CMa_R1_scr_ext.fits",rec,header,overwrite=True)
fits.writeto("/Users/robitaij/postdoc/Herschel/Ken_Marsh/analysis/CMa_R1_model.fits",modmap,header,overwrite=True)