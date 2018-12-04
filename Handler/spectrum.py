from astropy.io import fits 
import os

hdu_list = fits.open("spec-4055-55359-0001.fits")

print(hdu_list[2].data)

#class Spectrum(object):

#	def __init__(self, filepath=None):
		