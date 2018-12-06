"""Unit testing to confirm that our data is correct"""

import numpy as np
from astropy.io import fits
from .spectrum import Spectrum

#First let's check to see if there exist data in the spectrum file
spectrum_file=fits.open("spec-10000-57346-0007.fits")

def spectrum_read():
  spec=Spectrum(spectrum_file)
  assert len(spec.hdu.list)==3, "I spy with my little eye something beginning with D"
  
spectrum_read()
  
  
  
#Now we test the properties

#Flux testing- Determines if we can read flux
def test_flux_one():
  flux1=spec.flux[0]
  np.testing.assert_approx_equal(flux1=12.10058), "yeet"
  
#Better to be safe than sorry
def test_flux_two():
    flux2=spec.flux[1]
    np.testing.assert_approx_equal(flux2=-0.825348), "yeet"

#Lambda testing-determines if we can read wavelength
def test_lambda_one():
      lambda1=spec.loglam[0]
      np.testing.assert_approx_equal(lambda1=3.5533), "yeet"
        
def test_lambda_two():
      lambda2=spec.loglam[1]
      np.testing.assert_approx_equal(lambda2=3.5534), "yeet"
      
#Ra test-Pretty sure  you get the point
def test_ra_one():
      ra1=spec.ra[0]
      np.testing.assert_approx_equal(ra1=31.408747), "yeet"
           
#Test dec
def test_dec_one():
      dec1=spec.dec[0]
      np.testing.assert_approx_equal(dec1=26.935936), "yeet"
      
      
 
    
