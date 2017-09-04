from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
import reproject as proj
import numpy as np

#Read images
#----------------------------------

Qf = '/Users/jfrob/postdoc/GALFACTS/GALFACTS_S2_average_image_Q.fits'
Uf = '/Users/jfrob/postdoc/GALFACTS/GALFACTS_S2_average_image_U.fits'

Qhdu = fits.open(Qf)
Uhdu = fits.open(Uf)

#Modify headers
#----------------------------------

projhd = Qhdu[0].header
w = wcs.WCS(Qhdu[0].header)

projhd['CTYPE1'] = 'GLON-CAR'
projhd['CTYPE2'] = 'GLAT-CAR'

#Convert coordinates
val = w.all_pix2world([[projhd['CRPIX1'],projhd['CRPIX2'],0]],0)
Skyval = SkyCoord(ra=val[0,0]*u.degree, dec=val[0,1]*u.degree, frame='icrs')
Gval = Skyval.galactic

projhd['CRVAL1'] = Gval.l.degree
projhd['CRVAL2'] = Gval.b.degree

#Size of the projection
projhd['NAXIS1'] = 1074
projhd['NAXIS2'] = 1074

#Projection
#----------------------------------
Qproj = proj.reproject_exact(Qhdu,projhd)

fits.writeto('/Users/jfrob/postdoc/GALFACTS/GALFACTS_S2_average_image_Qgal.fits',Qproj,header=projhd)

projhd['OBJECT'] = Uhdu[0].header['OBJECT']
Uproj = proj.reproject_exact(Uhdu,projhd)
fits.writeto('/Users/jfrob/postdoc/GALFACTS/GALFACTS_S2_average_image_Ugal.fits',Uproj,header=projhd)