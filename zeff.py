import os
import numpy             as np
import pylab             as pl

import astropy.io.fits   as fits
import matplotlib.style  as style
import matplotlib.pyplot as plt


plt.figure(figsize=(15, 5))

root    = '/global/projecta/projectdirs/desi/datachallenge/redwood/'
truth   = root + '/targets/truth.fits'

dat     = fits.open(truth)[1]

tids    = dat.data['TARGETID'] 
tzs     = dat.data['TRUEZ']

print('\n')

for line, color, label in zip([1665.85, 2326.0, 2799.117, 3346.79, 3626., 3889.0, 4072.3, 4364.436], ['g', 'lime', 'chocolate', 'gold', 'c', 'r', 'm', 'k'], [r'OIII', 'CII', 'MgII', 'Ne V', 'OII', 'HeI', 'SII', 'OIII']):
  lo            = 4300. / line - 1.0
  hi            = 4500. / line - 1.0
  
  pl.fill_between(np.arange(lo, hi, 0.01), -2., 5., color=color, alpha=0.2)  
  plt.text((hi + lo) / 2., 4.7, label, {'ha': 'center'})
  
  print('{:.3f}  {:.3f}  {:.3f}'.format(line, lo, hi))


#  tiles   = ['10786', '8570', '8581', '8532', '8594', '8579', '8569', '8566', '8580']
tiles      = np.loadtxt('/global/homes/m/mjwilson/desi/bluebump/redwood_exposure_nums.txt', dtype=str)

#                                                                                                                                                                                                                               
for i, tile in enumerate(tiles): 
  redwood  = root + '/spectro/redux/redwood/spectra-64/{}/{}/zbest-64-{}.fits'.format(np.int(np.float(tile) / 100.), tile, tile)

  master   = os.environ['CSCRATCH'] + '/desi/bluebump/redwood/redrock/zbest-64-{}.fits'.format(tile)

  # S/N degraded according to throughput via IVAR. 
  degrade  = os.environ['CSCRATCH'] + '/desi/bluebump/degraded/redrock/zbest-64-{}.fits'.format(tile) 

  # Run in which actual throughput is degraded for quickspectra;  NOTE:  had no dynamic exposure times. 
  # dipthru  = os.environ['CSCRATCH'] + '/desi/bluebump/thrudip/redrock/zbest-64-{}.fits'.format(tile)
  
  ttargets = root + '/targets/{}/{}/truth-64-{}.fits'.format(np.int(np.float(tile) / 100.), tile, tile)
  
  result   = {} 

  for x, label in zip([redwood, master, degrade], ['Redwood', 'Master', 'Degraded', 'Dip']):
    result[label] = {}

    zbest         = fits.open(x)[1]
    ztids         = zbest.data['TARGETID']

    isin          = np.ones_like(ztids).astype(bool)
    
    '''
    if x == dipthru:
      ztids       = fits.open(ttargets)[1].data['TARGETID'][:5000]
      
      assert  np.all(ztids == np.sort(ztids))

      matched     = np.sort(result['Redwood']['ztids'][np.isin(result['Redwood']['ztids'], ztids)])
      isin        = np.isin(ztids, matched)

      result[label]['matched'] = np.isin(result['Redwood']['ztids'], matched)      
    '''

    #
    zz            = zbest.data['Z'][isin]
    zerr          = zbest.data['ZERR'][isin]
    zwarn         = zbest.data['ZWARN'][isin]
    ztids         = ztids[isin]

    # Match TIDS in spec. tile with target ids in truth table.
    index         = np.argsort(tids)
    sorted_tids   = tids[index]
    sorted_index  = np.searchsorted(sorted_tids, ztids)

    yindex        = np.take(index, sorted_index, mode="clip")
    mask          = tids[yindex] != ztids

    trusort       = np.ma.array(yindex, mask=mask)
    _             = tids[trusort]
    truzs         =  tzs[trusort]

    significance  = np.abs(tzs[trusort] - zz) / zerr
    significance  = np.log10(significance)

    result[label]['zbest']   = zbest
    result[label]['ztids']   = ztids
    result[label]['zz']      = zz
    result[label]['zerr']    = zerr
    result[label]['zwarn']   = zwarn
    result[label]['trusort'] = trusort
    result[label]['truzs']   = truzs
    result[label]['signif']  = significance

    print('\n\nHealpixel:  {}'.format(tile))                                                                                                                                                                                                                                                                   
    print('Number of targets: {}'.format(len(result[label]['ztids'])))
    print('Number of warnings for {}: {}'.format(label, np.count_nonzero(zwarn)))
    print('Number of redshifts in err by 1 sigma: {}'.format(np.count_nonzero(significance > 1.0)))
    print('Number of redshifts in err by 2 sigma: {}'.format(np.count_nonzero(significance > 2.0)))                                                                                                                                         print('Number of redshifts in err by 3 sigma: {}'.format(np.count_nonzero(significance > 3.0)))  

    print(ztids, _)
    
  #   
  # problem       = np.abs(result['Degraded']['zz'] - result['Master']['zz']) / result['Master']['zerr']
  # problem       = problem > 1.0

  # noproblem     = ~problem

  warning         = (result['Master']['zwarn'] > 0)
  
  pl.plot(tzs[result['Master']['trusort']][ warning], np.log10(np.abs(result['Degraded']['zz'] - result['Master']['zz']) / result['Master']['zerr'])[ warning], 'x', c='r', markersize=3)
  pl.plot(tzs[result['Master']['trusort']][~warning], np.log10(np.abs(result['Degraded']['zz'] - result['Master']['zz']) / result['Master']['zerr'])[~warning], 'x', c='k', markersize=3)
  
  '''
  problem         = np.abs(result['Dip']['zz'] - result['Master']['zz'][result['Dip']['matched']]) / result['Master']['zerr'][result['Dip']['matched']]
  problem         = problem > 1.0

  noproblem       = ~problem

  warning         = result['Dip']['zwarn'] == 0
  
  pl.plot(tzs[result['Master']['trusort']][result['Dip']['matched']][ warning], np.log10(np.abs(result['Dip']['zz'] - result['Master']['zz'][result['Dip']['matched']]) / result['Master']['zerr'][result['Dip']['matched']])[ warning], 'x', c='k', markersize=3)
  pl.plot(tzs[result['Master']['trusort']][result['Dip']['matched']][~warning], np.log10(np.abs(result['Dip']['zz'] - result['Master']['zz'][result['Dip']['matched']]) / result['Master']['zerr'][result['Dip']['matched']])[~warning], 'x', c='r', markersize=3)   
  '''
  
pl.xlim(-0.02, 1.8)
pl.ylim(-2.00, 5.0)
pl.xlabel(r'$z_{\rm{True}}$')
pl.ylabel(r"$\log_{10}(|z' - z| \ / \ z_{\rm{err}})$")
pl.legend(loc=1, frameon=True)

ax = pl.gca()

ax.grid(False)

plt.tight_layout()

pl.savefig('zeff.pdf')

print('\n\nDone.\n\n')
