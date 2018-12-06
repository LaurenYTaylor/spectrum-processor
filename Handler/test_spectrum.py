"""Unit testing to confirm that our data is correct"""

import numpy as np
from ..Spectrum import spectrum

#First let's check to see if there exist data in the spectrum file
spectrum_file="spec-3586-55181-0042.fits"

def spectrum_read():
  spec=Spectrum(spectrum_file)
  assert len(spec.hdu.list)==3, "I spy with my little eye something beginning with D"
  
  
#Now we test the properties

#Flux testing
def test_flux_one():
  #Tests to see if flux at [...] is True for x=....16.0137
  flux1=flux.spec[0]
  np.testing.assert_approx_equal(s.flux,16.0137)
  
#Better to be safe than sorry
def test_flux_two():
    flux2=flux.spec[1]
    np.testing.assert_approx_equal(s.flux,xxx)

 #Lambda testing
def test_lambda_one():
      lambda1=loglam.spec[0]
      np.testing.assert_approx_equal(s.loglam,xxx)
        
def test_lambda_two():
      lambda2=loglam.spec[1]
      np.testing.assert_approx_equal(s.loglam,xxx)
      

   
    
