import numpy           as     np
import pylab           as     pl
import astropy.io.fits as     fits

from   scipy           import interpolate
from   astropy.table   import Table


file  = open('redwood_spectra_files.txt', mode = 'r')
files = file.readlines()
files = [x.replace('\n', '') for x in files]
file.close()

wave  = fits.open('data/throughput/thru-b-dip.fits')[1].data['wavelength']
_thru = fits.open('data/throughput/thru-b-dip.fits')[1].data['throughput']

thru  = interpolate.interp1d(wave, _thru, fill_value=1.0, bounds_error=False)

for x in files:
  _      = x.split('/')[-1]

  spec   = fits.open(x + '/spectra-64-{}.fits'.format(_))
  f      = thru(spec['B_WAVELENGTH'].data) 

  spec['B_IVAR'].data = f * f * spec['B_IVAR'].data
  
  spec.writeto('/global/cscratch1/sd/mjwilson/desi/bluebump/degraded/spectra-64-{}.fits'.format(_), overwrite=True)
