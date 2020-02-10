import pylab as pl
import numpy as np
import astropy.io.fits as fits

from astropy.table import Table


root = '/global/projecta/projectdirs/desi/datachallenge/redwood/targets/103/10323/'

wave = fits.open(root + 'truth-64-10323.fits')[2].data  
flux = fits.open(root + 'truth-64-10323.fits')[3].data

for x in flux[:10]: 
  pl.plot(wave, x)

pl.savefig('redwood_tflxs.pdf')
