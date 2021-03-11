# PYWAVAN

## Wavelet transforms functions for python

The main function of this repository is **fan_trans()**, which can be found in *wavan.py*

The function **fan_trans()** performs the Fan wavelet transform as described by [Kirby 2005](http://adsabs.harvard.edu/abs/2005CG.....31..846K).

Additional options allow the user to perform a Multiscale non-Gaussian segmentation as descirbed by [Robitaille et al. 2014](https://ui.adsabs.harvard.edu/abs/2014MNRAS.440.2726R/abstract) and [Robitaille et al. 2019](https://ui.adsabs.harvard.edu/abs/2019arXiv190511492R/abstract).

Here is a non-hexaustive description for some other functions:

* The function **powspec()** calculates the classical Fourier Power Spectrum of an image where the Fourier coefficients are averaged over the azimuthal angles as a function of the wavenumber.

* The function **gauss_beam()** allows one to model the transfer function of a telescope by a Gaussian function.

* Functions **subfits()**, **congrid()** (in *rebin.py*) and **imsmooth()** allow respectively to cut a fits file, rebin a map and smooth a map.

* Finally, the function **fbm2d()** in *noisegen.py* allows one to create a fractional Brownian motion map, i.e. a fractal image, respecting a specific power law.

**powspec()**, **fbm2d()** and **apodize()** functions are inspired by the [mamdlib](https://www.ias.u-psud.fr/pperso/mmiville/mamdlib.html) IDL library.

## Installation

**Dependencies**

`pip install tqdm`

`pip install astropy`

**pywavan installation**

Clone this repository:

`git clone https://github.com/jfrob27/pywavan.git`

and in the *pywavan* directory

`pip install .`

or directly

`pip install git+https://github.com/jfrob27/pywavan.git`

## Tutorials

An introduction Notebook tutorial is available [here](https://nbviewer.jupyter.org/github/jfrob27/pywavan/blob/master/tutorial/tuto_pywavan.ipynb) and in the *tutorial* folder.

Pywavan is also demonstrated on real astrophysical data along some [Turbustat](https://github.com/Astroua/TurbuStat) functions in the following [Notebook](https://nbviewer.jupyter.org/github/jfrob27/workshop/blob/main/multi-scale_analysis.ipynb). This Notebook is also designed to work directly in [Google Colab](https://colab.research.google.com/github/jfrob27/workshop/blob/main/multi-scale_analysis.ipynb).