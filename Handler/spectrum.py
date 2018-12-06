from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy import constants as const
import os
from matplotlib import pyplot as plt
import numpy as np

class Redshift(object):

	def __init__(self, z):
		self.H_0 = 67.77*((u.km/u.s)/u.Mpc)
		self.z=z
	
	@property
	def velocity(self):
		velocity = self.z*const.c
		return velocity
	
	@property
	def distance(self):
		distance = self.velocity/self.H_0
		distance = distance.decompose()
		distance = distance.to(u.kpc)
		return distance


class Spectrum(object):
	def __init__(self, filepath):
		self.filepath = filepath
	
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
		
	@property 
	def min_lambda(self):
		t = Table.read(self.filepath, hdu=2)
		min = t['WAVEMIN'].data[0]
		return min
	
	@property 
	def max_lambda(self):
		t = Table.read(self.filepath, hdu=2)
		max = t['WAVEMAX'].data[0]
		return max
	
	@property
	def redshift(self):
		t = Table.read(self.filepath, hdu=2)
		redshift = Redshift(t['Z'].data[0])
		return redshift
		
	def display_headers(self, header_num):
		t = Table.read(self.filepath, hdu=header_num)
		print(t)
		
	def display_info(self):
		hdu = fits.open(self.filepath)
		return hdu.info()
		
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
print(spec.redshift.distance)


