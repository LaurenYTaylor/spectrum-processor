"""Unit testing to confirm that our data is correct"""

import numpy as np
from astropy.io import fits
from Spectrum import spectrum

#First let's check to see if there exist data in the spectrum file
spectrum_file=fits.open("spec-3586-55181-0042.fits")

def spectrum_read():
  spec=Spectrum(spectrum_file)
  assert len(spec.hdu.list)==3, "I spy with my little eye something beginning with D"
  
  
#Now we test the properties

#Flux testing- Determines if we can read flux
def test_flux_one():
  flux1=spec.flux[0]
  np.testing.assert_approx_equal(flux1=16.0137)
  
#Better to be safe than sorry
def test_flux_two():
    flux2=spec.flux[1]
    np.testing.assert_approx_equal(flux2=xxx)

#Lambda testing-determines if we can read wavelength
def test_lambda_one():
      lambda1=spec.loglam[0]
      np.testing.assert_approx_equal(lambda1=xxx)
        
def test_lambda_two():
      lambda2=spec.loglam[1]
      np.testing.assert_approx_equal(lambda2=xxx)
      
#Ra test-Pretty sure  you get the point
def test_ra_one():
      ra1=spec.ra[0]
      np.testing.assert_approx_equal(ra1,xxx)
      
def test_ra_two():
      ra2=spec.ra[1]
      np.testing.assert_approx_equal(ra2=xxx)
      
#Test dec
def test_dec_one():
      dec1=spec.dec[0]
      np.testing.assert_approx_equal(dec1=xxx)
      
def test_dec_two():
      dec2=spec.dec[1]
      np.testing.assert_approx_equal(dec2=xxx)
      
 
    
