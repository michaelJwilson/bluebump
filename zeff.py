import os
import numpy as np
import pylab as pl

import astropy.io.fits as fits


root  = '/global/projecta/projectdirs/desi/datachallenge/redwood/'
truth = root + '/targets/truth.fits'

oput  = os.environ['CSCRATCH'] + '/desi/bluebump/zbest-64-10323.fits'

dat   = fits.open(truth)[1]

print(dat.data)

tids  = dat.data['TARGETID'] 
tzs   = dat.data['TRUEZ']

zbest = fits.open(oput)[1]
ztids = zbest.data['TARGETID']
zz    = zbest.data['Z']
zerr  = zbest.data['ZERR']
zwarn = zbest.data['ZWARN']

index         = np.argsort(tids)
sorted_tids   = tids[index]
sorted_index  = np.searchsorted(sorted_tids, ztids)

yindex        = np.take(index, sorted_index, mode="clip")
mask          = tids[yindex] != ztids

result        = np.ma.array(yindex, mask=mask)

pl.plot(tzs[result], zz, '.')

pl.savefig('zeff.pdf')


