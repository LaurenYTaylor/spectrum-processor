"""Unit testing to confirm that our data is correct"""

import math 
from spectrum import spectrum

#First let's check to see if there exist data in the spectrum file
spectrum_file="sdss.spec-xxxx-xxxxx-xxxx.fits"

def test_read_flux():
  flux=Spectrum("spectrum_file")
  assert len(flux.hdu.list)>4, "I spy with my little eye something beginning with D"
  
def test_read_lambda():
  lambda=Spectrum("spectrum_file")
  assert len(lambda.hdu.list)>4 "The D stands for data"
  
#Now we test if specific values correlate with what the data is

def test_flux_one():
  #Tests to see if flux at [...] is True for x=....
  flux1=spectrum(#insert flux row number here)
  assert flux1 == #insert value at flux row, "Data lines up with input"
    
#Better to be safe than sorry
def test_flux_two():
    flux2=spectrum(#insert flux row number here)
    assert flux2 == #insert value at row number, "Data lines up with input"
      
def test_lambda_one():
      lambda1=spectrum(#insert loglam row number here)
      assert lambda1 == #insert value at row number, "Data lines up with input"
        
def test_lambda_two():
      lambda2=spectrum(#insert loglam row number here)
      assert lambda2 == #insert value at row number, "Data lines up with input"
       
   
    
