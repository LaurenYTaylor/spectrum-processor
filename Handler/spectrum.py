from astropy.io import fits
from astropy.table import Table
import os
from matplotlib import pyplot as plt
import numpy as np

# hdu_info = fits.open("spec-10000-57346-0007.fits")
# t1 = Table.read("spec-10000-57346-0007.fits", hdu=1)
# t2 = Table.read("spec-10000-57346-0007.fits", hdu=2)

# print(t2['PLUG_RA'].data)
# print(t2['PLUG_DEC'].data)
# print(t2['OBJTYPE'].data)
# print(t2['CLASS'].data)
# print(t2['SUBCLASS'].data)


class Spectrum(object):
	def __init__(self, filepath):
		self.filepath = filepath
		#self.targeted = None
		#self.obs_date = None
		#self.z = None
		#self.distance = None
		#self.max_flux = None
		#self.obj_class = None
		#self.obj_subclass = None
	
	@property
	def flux(self):
		t = Table.read(self.filepath, hdu=1)
		flux = t['flux'].data
		return flux
			
	@property
	def loglam(self):
		t = Table.read(self.filepath, hdu=1)
		loglam = t['loglam'].data
		return loglam	
	
	@property
	def ra(self):
		t = Table.read(self.filepath, hdu=2)
		ra = t['PLUG_RA'].data[0]
		return ra
	
	@property
	def dec(self):
		t = Table.read(self.filepath, hdu=2)
		dec = t['PLUG_DEC'].data[0]
		return dec
	
	def display_headers(self, header_num):
		t = Table.read(self.filepath, hdu=header_num)
		print(t)
		
	def display_info(self):
		hdu = fits.open(self.filepath)
		print(hdu.info())
		
	
spec = Spectrum("spec-10000-57346-0007.fits")
print(spec.ra)
	