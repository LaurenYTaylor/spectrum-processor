from astropy.io import fits
from astropy.table import Table
import os
from matplotlib import pyplot as plt
import numpy as np

# hdu_info = fits.open("spec-10000-57346-0007.fits")
# t1 = Table.read("spec-10000-57346-0007.fits", hdu=1)
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
		
	def plot_spectrum(self, show):
		plt.rc('text', usetex=True)
		plt.plot(self.loglam, self.flux)
		plt.xlabel(r"Log $\lambda$") 
		plt.ylabel(r"Flux $\frac{erg.\AA}{s.cm^2}$")
		if not os.path.isdir("Plots"):
			os.mkdir('Plots')
		plt.savefig("Plots/"+self.filepath[:-5]+".png")
		if(show==1 or show==True):
			plt.show()
		
	
spec = Spectrum("spec-10000-57346-0007.fits")
#print(spec.ra)
#spec.display_headers(1)
#spec.display_info()	
spec.plot_spectrum(show=1)
print(np.sum(spec.flux))