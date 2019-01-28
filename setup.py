from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='pywavan',
      version='0.1',
      description='Functions to perform the Multiscale non-Gaussian segmentation (MnGSeg) and essential tools to perform a power spectrum analysis',
      url='https://github.com/jfrob27/pywavan',
      author='J.-F. Robitaille',
      author_email='jean-francois.robitaille@univ-grenoble-alpes.fr',
	  long_description=long_description,
	  long_description_content_type="text/markdown",
      license='New BSD',
      packages=['pywavan'],
      zip_safe=False)