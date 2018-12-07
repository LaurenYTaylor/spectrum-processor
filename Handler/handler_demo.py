import argparse
from spectrum import Spectrum
parser = argparse.ArgumentParser(description="A script to demonstrate an SDSS spectrum fits file handler.",\
								usage="handler_demo.py --filename spec-xxxx-xxxxx-xxxx.fits")

parser.add_argument("-f", "--filename", help="The spectrum fits file to read.", default=10)
args = parser.parse_args()

spectrum_file = args.filename

spectrum = Spectrum(spectrum_file)

redshift = spectrum.redshift.z
velocity = spectrum.redshift.velocity
distance = spectrum.redshift.distance
print(f"Redshift: {redshift}")
print(f"Velocity (from redshift): {velocity}")
print(f"Distance (from redshift): {distance}")

print(spectrum.luminosity)
