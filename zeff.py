import os
import numpy             as np
import pylab             as pl

import astropy.io.fits   as fits
import matplotlib.style  as style
import matplotlib.pyplot as plt


# style.use('bmh')

plt.figure(figsize=(15, 5))

root    = '/global/projecta/projectdirs/desi/datachallenge/redwood/'
truth   = root + '/targets/truth.fits'

dat     = fits.open(truth)[1]

tids    = dat.data['TARGETID'] 
tzs     = dat.data['TRUEZ']

print('\n')

for line, color, label in zip([1216., 2799.117, 3626., 4000.], ['g', 'k', 'c', 'm'], [r'Ly-$\alpha$', 'Mg-II', 'OII', r'H-$\alpha$']):
  lo            = 4300. / line - 1.0
  hi            = 4500. / line - 1.0
  
  pl.fill_between(np.arange(lo, hi, 0.01), -2., 5., color=color, alpha=0.2)  
  plt.text((hi + lo) / 2., 4.7, label, {'ha': 'center'})
  
  print('{:.3f}  {:.3f}  {:.3f}'.format(line, lo, hi))

#                                                                                                                                                                                                                               
for i, tile in enumerate(['10786', '8570', '8581', '8532', '8594', '8579', '8569', '8566', '8580']): 
  degrade = os.environ['CSCRATCH'] + '/desi/bluebump/degraded/redrock/zbest-64-{}.fits'.format(tile)
  redwood = os.environ['CSCRATCH'] + '/desi/bluebump/redwood/redrock/zbest-64-{}.fits'.format(tile)
  
  for x, color, label in zip([redwood, degrade], ['r', 'b'], ['Original', 'Degraded']):
    zbest         = fits.open(x)[1]
    ztids         = zbest.data['TARGETID']
    zz            = zbest.data['Z']
    zerr          = zbest.data['ZERR']
    zwarn         = zbest.data['ZWARN']

    index         = np.argsort(tids)
    sorted_tids   = tids[index]
    sorted_index  = np.searchsorted(sorted_tids, ztids)

    yindex        = np.take(index, sorted_index, mode="clip")
    mask          = tids[yindex] != ztids

    result        = np.ma.array(yindex, mask=mask)

    significance  = np.abs(tzs[result] - zz) / zerr
    significance  = np.log10(significance)
  
    if color == 'r':
      shifts      = np.random.normal(loc=0.0, scale=0.05, size=len(significance))

    if i > 0:
      label = ''
      
    pl.plot(tzs[result] + shifts, significance + shifts, 'x', c=color, markersize=3, label=label)

    print('\n\nTile:  {}'.format(tile))
    print('Number of warnings for {}: {}'.format(label, np.count_nonzero(zwarn)))
    print('Number of redshifts in err by 1 sigma: {}'.format(np.count_nonzero(significance > 1.0)))
    print('Number of redshifts in err by 2 sigma: {}'.format(np.count_nonzero(significance > 2.0)))
    print('Number of redshifts in err by 3 sigma: {}'.format(np.count_nonzero(significance > 3.0)))
  
pl.xlim(0.0, 1.8)
pl.ylim(-2., 5.0)
pl.xlabel(r'$z_{\rm{True}}$')
pl.ylabel(r'$\log_{10}(|\hat z - z_{\rm{True}}| \ / \ z_{\rm{err}})$')
pl.legend(loc=1, frameon=True)

ax = pl.gca()

ax.grid(False)

plt.tight_layout()

pl.savefig('zeff.pdf')

print('\n\nDone.\n\n')
