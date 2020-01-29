import os
import numpy as np
import pylab as pl

import astropy.io.fits   as fits
import matplotlib.style  as style
import matplotlib.pyplot as plt


style.use('bmh')

root    = '/global/projecta/projectdirs/desi/datachallenge/redwood/'
truth   = root + '/targets/truth.fits'

dat     = fits.open(truth)[1]

tids    = dat.data['TARGETID'] 
tzs     = dat.data['TRUEZ']

#
degrade = os.environ['CSCRATCH'] + '/desi/bluebump/degraded/redrock/zbest-64-8579.fits'
redwood = os.environ['CSCRATCH'] + '/desi/bluebump/redwood/redrock/zbest-64-8579.fits'

shift   = {'r': 0.025, 'b': 0.0}

for x, color, label in zip([degrade, redwood], ['r', 'b'], ['Degraded', 'Original']):
  zbest = fits.open(x)[1]
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

  pl.plot(tzs[result] + shift[color], zz, 'x', c=color, markersize=3, label=label)

pl.xlabel(r'$z_{\rm{True}}$')
pl.ylabel(r'$\hat z$')
pl.legend(loc=2, frameon=False)

ax = pl.gca()

ax.grid(False)

plt.tight_layout()

pl.savefig('zeff.pdf')


