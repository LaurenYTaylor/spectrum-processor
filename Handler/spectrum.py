from astropy.io import fits
from astropy.table import Table, Column
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



class Geometry(object):
	
	def __init__(self, r):
		
		self.r = r

	def area(self):

		area = 4*np.pi*self.r*self.r
		return area

class Luminosity(object):

	def __init__(self, flux, distance):
		self.flux=flux
		self.distance=distance
	
	@property
	def luminosity(self):
		g = Geometry(self.distance)
		luminos = self.flux*g.area()
		luminos = luminos.decompose()
		return luminos

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
	
	@property
	def luminosity(self):
		lum = Luminosity(self.flux, self.redshift.distance)
		return lum.luminosity
	
	@property 
	def object_type(self):
		t = Table.read(self.filepath, hdu=2)
		object_type = t['CLASS'].data[0]
		return object_type
	
	def display_headers(self, header_num):
		t = Table.read(self.filepath, hdu=header_num)
		print(t)
		
	def display_info(self):
		hdu = fits.open(self.filepath)
		return hdu.info()
		
	def plot_spectrum(self, show=True, plotlines='all'):

		fig, ax = plt.subplots(1)

		if plotlines == None:
			pass
		elif plotlines == 'all':
			lines = SpectralLines(self.filepath)
			lines.plot_all_lines(ax)
			w, n = lines.get_all_lines()
			testw = np.log10(w)
			print("\n---------------\nPlotting lines:\n---------------")
			for i in range(len(w)):
				if np.log10(self.min_lambda) < testw[i] < np.log10(self.max_lambda):
					print(f"{n[i]} {w[i]:.2f} angstroms: {np.log10(w[i]):.2f} log lambda")
		else:
			lines = SpectralLines(self.filepath)
			lines.plot_some_lines(ax, plotlines)
			w, n = lines.get_lines(plotlines)
			testw = np.log10(w)
			print("\n---------------\nPlotting lines:\n---------------")
			for i in range(len(w)):
				if np.log10(self.min_lambda) < testw[i] < np.log10(self.max_lambda):
					print(f"{n[i]} {w[i]:.2f} angstroms: {np.log10(w[i]):.2f} long lambda")

		plt.rc('text', usetex=True)
		plt.plot(self.loglam, self.flux, 'k-', lw=0.5)
		plt.xlim(np.log10(self.min_lambda), np.log10(self.max_lambda))
		plt.xlabel(r"Log $\lambda$") 
		plt.ylabel(r"Flux $\frac{erg.\AA}{s.cm^2}$")
		plt.title("Flux vs. Log Wavelength")
		if not os.path.isdir("Plots"):
			os.mkdir('Plots')
		plt.savefig("Plots/"+self.filepath[:-5]+".png")
		if(show==1 or show==True):
			plt.show()


class SpectralLines(object):

	def __init__(self, filepath):
		self.filepath = filepath

	@property
	def lines(self):
		t = Table.read(self.filepath, hdu=3)
		return t['LINEWAVE'].data

	@property
	def linenames(self):
		t = Table.read(self.filepath, hdu=3)
		names = []
		for i, val in enumerate(t['LINENAME'].data):
			names.append(val.decode("utf-8").strip('[] 0123456789'))
		return names

	@property
	def colours(self):
		c = ['#2B0642', '#2B0642','#2B0642','#2B0642','#2B0642','#2B0642','#2B0642','#2B0642', '#740080', '#780088', '#8100a9', '#7e00db', '#2800ff', '#1d00ff', '#00a5ff', '#00efff', '#00ffc0', '#00ff87', '#85ff00', '#bdff00', '#f3ff00', '#ffe600', '#ff4f00', '#ff3400', '#ff4b00', '#ff3000', '#ff0000', '#ff0000', '#ff0000', '#ff0000', '#ff0000', '#e60000']
		return c

	def to_log(self, lam):
		return np.log10(lam)
		
	def display_info(self):
		hdu = fits.open(self.filepath)
		print(hdu.info())

	def get_one_line(self, lam):
		wavelength = self.lines[np.where(np.isclose(lam, self.lines, atol=5))][0]
		name = self.linenames[np.where(np.isclose(lam, self.lines, atol=5))[0][0]]
		return wavelength, name

	def get_lines(self, lams):
		if isinstance(lams, int):	
			lams = [lams]
		wavelengths = []
		names = []
		for i in range(len(lams)):
			try:
				temp_w, temp_n = self.get_one_line(lams[i])
				wavelengths.append(temp_w)
				names.append(temp_n)
			except IndexError:
				a=0
				#print(f"No spectral line on file within five angstroms either side of {lams[i]} angstroms")
		return wavelengths, names

	def get_all_lines(self):
		return self.lines, self.linenames

	def plot_some_lines(self, ax, lams):
		if isinstance(lams, int):
			lams = [lams]
		w, n = self.get_lines(np.sort(lams))
		lines = self.to_log(w)
		c = np.array(self.colours)[np.where(np.isin(self.lines, w))]
		plt.vlines(lines, 0, 1, transform=ax.get_xaxis_transform(), colors=c, linestyle='--', alpha=0.75)

	def plot_all_lines(self, ax):
		lines = self.to_log(self.lines)
		plt.vlines(lines, 0, 1, transform=ax.get_xaxis_transform(), colors=self.colours, linestyle='--')