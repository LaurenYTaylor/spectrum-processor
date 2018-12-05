"""Unit testing to confirm that our data is correct"""

import math 
from spectrum import spectrum

#First let's check to see if there exist data in the spectrum file
spectrum_file="spec-xxxx-xxxxx-xxxx.fits"

def test_read_flux():
  flux=Spectrum("spectrum_file")
  assert len(flux.hdu.list)>4, "I spy with my little eye something beginning with D"
  
def test_read_lambda():
  lambda=Spectrum("spectrum_file")
  assert len(lambda.hdu.list)>4 "The D stands for data"
  
  
  

def test_spectrumone():
  #Tests to see if flux at [...] is True for x=....
  spec1=spectrum(...)
  assert
