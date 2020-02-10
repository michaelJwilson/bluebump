import numpy           as     np
import pylab           as     pl
import astropy.io.fits as     fits

from   scipy           import interpolate
from   astropy.table   import Table


file  = open('redwood_spectra_files.txt', mode = 'r')
files = file.readlines()
files = [x.replace('\n', '') for x in files]
file.close()

wave, _thru = np.loadtxt('thru-dip.txt', unpack=True)
thru        = interpolate.interp1d(wave, _thru, fill_value=1.0, bounds_error=False)

for x in files:
  _    = x.split('/')[-1]

  spec = fits.open(x + '/spectra-64-{}.fits'.format(_))
  f    = thru(spec['B_WAVELENGTH'].data) 

  # pl.plot(spec['B_WAVELENGTH'].data,  spec['B_IVAR'].data[0,:], 'b')
   
  spec['B_IVAR'].data = f * f * spec['B_IVAR'].data
  
  spec.writeto('/global/cscratch1/sd/mjwilson/desi/bluebump/degraded/spectra-64-{}.fits'.format(_), overwrite=True)

  # Check.
  # dspec = fits.open('/global/cscratch1/sd/mjwilson/desi/bluebump/degraded/spectra-64-{}.fits'.format(_))
  # pl.plot(spec['B_WAVELENGTH'].data, dspec['B_IVAR'].data[0,:], 'r')
  # pl.savefig('divar.pdf')
  # break
