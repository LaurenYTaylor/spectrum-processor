import argparse
from spectrum import Spectrum
parser = argparse.ArgumentParser(description="A script to demonstrate an SDSS spectrum fits file handler.",\
								usage="handler_demo.py --filename spec-xxxx-xxxxx-xxxx.fits")

parser.add_argument("-f", "--filename", help="The spectrum fits file to read.", default=10, required=True)
args = parser.parse_args()

spectrum_file = args.filename

spectrum = Spectrum(spectrum_file)

redshift = spectrum.redshift.z
velocity = spectrum.redshift.velocity
distance = spectrum.redshift.distance
print(f"\n--------\nRedshift\n--------\n{redshift}")
print(f"\n------------------------\nVelocity (from redshift)\n------------------------\n{velocity}")
print(f"\n------------------------\nDistance (from redshift)\n------------------------\n{distance}")
print("\n-----------\nLuminosity\n-----------")
for l in spectrum.luminosity[:5]:
	print(l)
#print("\n----\nInfo\n----")
#print(spectrum.display_info())

#print("\n---------\nHeader 1\n---------")
#print(spectrum.display_headers(1))

#spectrum.plot_spectrum(show=True, plotlines=None) # for no lines
#spectrum.plot_spectrum(show=True, plotlines=[5578, 5877]) # the two lines in file "spec-10000-57346-0007.fits"
spectrum.plot_spectrum() # alternative to show default line plot; can also show plotlines=None
