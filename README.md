# PYWAVAN

## Wavelet transforms functions for python

The main function of this repository is **fan_trans()**, which can be found in *wavan.py*

The function **fan_trans()** performs the Fan wavelet transform as described by [Kirby 2005](http://adsabs.harvard.edu/abs/2005CG.....31..846K).

Additional options allow the user to perform a Multiscale non-Gaussian segmentation as descirbed by [Robitaille et al. 2014](http://adsabs.harvard.edu/abs/2014MNRAS.440.2726R) and Robitaille et al. 2019 (in preparation).

Here is a non-hexaustive description for some other functions:

* The function **powspec()** calculates the classical Fourier Power Spectrum of an image where the Fourier coefficients are averaged over the azimuthal angles as a function of the wavenumber.

* The function **gauss_beam()** allows one to model the transfer function of a telescope by a Gaussian function.

* Functions **subfits()**, **congrid()** (in *rebin.py*) and **imsmooth()** allow respectively to cut a fits file, rebin a map and smooth a map.

* Finally, the function **fbm2d()** in *noisegen.py* allows one to create a fractional Brownian motion map, i.e. a fractal image, respecting a specific power law.

Many of these functions are strongly inspired by the [mamdlib](https://www.ias.u-psud.fr/pperso/mmiville/mamdlib.html) IDL library.

## Installation

Clone this repository:

`git clone https://github.com/jfrob27/pywavan.git`

and add the folder path to your `PYTHONPATH`.

OR

`pip install -i https://test.pypi.org/simple/ pywavan`

(however pywavan is still at is early stage of delivery/testing)