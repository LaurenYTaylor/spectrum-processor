from astropy.io import fits
from astropy.table import Table
import os
from matplotlib import pyplot as plt
import numpy as np

hdu_info = fits.open("spec-10000-57346-0007.fits")
t1 = Table.read("spec-10000-57346-0007.fits", hdu=1)
t2 = Table.read("spec-10000-57346-0007.fits", hdu=2)

print(t2['PLUG_RA'].data)
print(t2['PLUG_DEC'].data)
print(t2['OBJTYPE'].data)
print(t2['CLASS'].data)
print(t2['SUBCLASS'].data)


class Spectrum(object):
	def __init__(self, filepath):
		filepath = filepath
		self.flux = None
		self.loglam = None
		self.ra = None
		self.dec = None
		self.targeted = None
		self.obs_date = None
		self.z = None
		self.distance = None
		self.max_flux = None
		self.obj_class = None
		self.obj_subclass = None
	
	@getter
	def flux(self):
		t = Table.read(self.filepath, hdu=1)
		flux = t['flux'].data
		return flux
	
	
	def display_headers(self, header_num):
		t = Table.read(self.filepath, hdu=header_num)
		print(t)
		
	def display_info(self):
		hdu = fits.open(self.filepath)
		print(hdu.info())
		
	
		