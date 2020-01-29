import numpy as np
import pylab as pl
import astropy.io.fits as fits

from   scipy         import interpolate
from   astropy.table import Table


dat     = np.loadtxt('thru-dip.txt')

pl.plot(dat[:,0], dat[:,1])

##  Blue limit of r camera.
pl.axvline(x=5564.000, ymin=0, ymax=1, c='k')

print('Relative throughput: {} to {}'.format(dat[:,1].min(), dat[:,1].max()))

affected = (dat[:,1] < 1.00) | (dat[:,1] > 1.00)
awave    =  dat[:,0][affected]

print('Affected range: {} to {}'.format(awave.min(), awave.max()))

dthru    = interpolate.interp1d(dat[:,0], dat[:,1], fill_value=1.0, bounds_error=False)


##  DESI b throughput
dat      = Table(fits.open('data/throughput/thru-b.fits')[1].data)

pl.plot(dat['wavelength'], dthru(dat['wavelength']), 'r')


##  New throughput
dat['throughput'] *= dthru(dat['wavelength'])

dat.write('data/throughput/thru-b-dip.fits', format='fits', overwrite=True)

pl.savefig('dip.pdf')

print('\n\nDone.\n\n')
