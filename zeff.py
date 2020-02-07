import os
import pickle
import numpy             as np
import pylab             as pl
import pandas            as pd

import astropy.io.fits   as fits
import matplotlib.style  as style
import matplotlib.pyplot as plt

from   astropy.table     import Table


plt.figure(figsize=(15, 5))

def significance(truzs, zz, zerr):
  significance = np.abs(truzs - zz) / zerr
  significance = np.log10(significance)

  return  significance


root    = '/global/projecta/projectdirs/desi/datachallenge/redwood/'
truth   = root + '/targets/truth.fits'

truth   = Table.read(truth, format='fits')
# truth = truth.to_pandas() 

tids    = truth['TARGETID'] # .to_numpy()
tzs     = truth['TRUEZ']    # .to_numpy()

del truth['MOCKID']
del truth['SEED']
del truth['TEFF']
del truth['LOGG']
del truth['FEH']

# Nanomaggies to mags.
for x in ['FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2']:
  ff       = truth[x]  
  truth[x] = 22.5 - 2.5 * np.log10(ff)
  
print(len(truth))
print(truth)

print('\n')

for line, color, label in zip([1549.5, 1908.7, 2799.117, 3727., 3889.0, 4072.3, 4364.436], ['k', 'g', 'chocolate', 'c', 'r', 'm'], ['CIV', 'CIII', 'MgII', 'OII', 'HeI', 'SII', 'OIII']):
  lo            = 4300. / line - 1.0
  hi            = 4500. / line - 1.0
  
  pl.fill_between(np.arange(lo, hi, 0.01), -2., 5., color=color, alpha=0.2)  
  plt.text((hi + lo) / 2., 4.7, label, {'ha': 'center'})
  
  print('{:.3f}  {:.3f}  {:.3f}'.format(line, lo, hi))



tiles      = np.loadtxt('/global/homes/m/mjwilson/desi/bluebump/redwood_exposure_nums.txt', dtype=str)
labels     = ['BAD Z', 'LRG', 'ELG', 'QSO'] 

results    = []

for i, tile in enumerate(tiles): 
  # Original redwood run files.
  redwood  = root + '/spectro/redux/redwood/spectra-64/{}/{}/zbest-64-{}.fits'.format(np.int(np.float(tile) / 100.), tile, tile)

  # Rerun of redwood spectra with redrock master: 02/07/20.
  master   = os.environ['CSCRATCH'] + '/desi/bluebump/redwood/redrock/zbest-64-{}.fits'.format(tile)

  # S/N degraded according to throughput via IVAR; redrock master.  
  degrade  = os.environ['CSCRATCH'] + '/desi/bluebump/degraded/redrock/zbest-64-{}.fits'.format(tile) 

  # Run in which actual throughput is degraded for quickspectra;  NOTE:  had no dynamic exposure times so not directly comparable to redwood. 
  # dipthru = os.environ['CSCRATCH'] + '/desi/bluebump/thrudip/redrock/zbest-64-{}.fits'.format(tile)

  # Corresponding Truth file.
  ttargets = root + '/targets/{}/{}/truth-64-{}.fits'.format(np.int(np.float(tile) / 100.), tile, tile)
  
  result   = {} 
  
  # 'Dip'
  for x, label in zip([redwood, master, degrade], ['Redwood', 'Master', 'Degraded']):
    # Load a spectroscopic tile.
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

    # TARGETIDS is the spec. tile.
    ztids         = ztids[isin]

    # Match TIDS in spec. tile with target ids in truth table with a sorted search.
    index         = np.argsort(tids)
    sorted_tids   = tids[index]
    sorted_index  = np.searchsorted(sorted_tids, ztids)

    yindex        = np.take(index, sorted_index, mode="clip")
    mask          = tids[yindex] != ztids

    trusort       = np.ma.array(yindex, mask=mask)

    # Spec. - matched TARGETIDs in truth. 
    _             = tids[trusort][~trusort.mask]

    # Corresponding redshifts.
    truzs         =  tzs[trusort][~trusort.mask]

    # Objects that weren't in the truth tables.  Standards?
    zzmask        = np.isin(ztids, _)

    # Check that TID matches between truth table and spec. tile was successful.
    assert  np.all(ztids[zzmask] == _)

    # Save to a dict.
    result[label]['zbest']   = zbest
    result[label]['ztids']   = ztids
    result[label]['zz']      = zz
    result[label]['zerr']    = zerr
    result[label]['zwarn']   = zwarn
    result[label]['zzmask']  = zzmask
    result[label]['trusort'] = trusort
    result[label]['truzs']   = truzs
    result[label]['signif']  = significance(truzs, zz[zzmask], zerr[zzmask])

    print('\n\nHealpixel:  {}'.format(tile))
    print('Number of targets: {}'.format(len(result[label]['ztids'])))
    print('Number of warnings for {}: {}'.format(label, np.count_nonzero(zwarn)))
    print('Number of redshifts in err by 1 sigma: {}'.format(np.count_nonzero(result[label]['signif'] > 1.0)))
    print('Number of redshifts in err by 2 sigma: {}'.format(np.count_nonzero(result[label]['signif'] > 2.0)))
    print('Number of redshifts in err by 3 sigma: {}'.format(np.count_nonzero(result[label]['signif'] > 3.0)))
    
  # Check that no funny sorting happened between the Master and Degraded runs. 
  assert  np.all(result['Degraded']['ztids'] == result['Master']['ztids'])

  #
  zzmask       =  result['Master']['zzmask']
  warning      = (result['Master']['zwarn'][zzmask] > 0)

  sig          =  significance(result['Degraded']['zz'], result['Master']['zz'], result['Master']['zerr'])[zzmask][~warning]
  band         =  sig > 0.0

  mids         =  result['Master']['ztids'][zzmask][~warning]

  #  
  isin         =  np.isin(truth['TARGETID'], mids)
  _            =  Table(truth[isin], copy=True)

  _.sort('TARGETID')

  inds         =  np.argsort(mids)

  _['SIG']     =  sig[inds]

  select       = (_['SIG'] > 0.0)

  print()
  print(_[select])

  is_elg       = [x.strip() == 'ELG' for x in _['TEMPLATETYPE']]
  is_lrg       = [x.strip() == 'LRG' for x in _['TEMPLATETYPE']]
  is_qso       = [x.strip() == 'QSO' for x in _['TEMPLATETYPE']]

  try:
     pl.plot(truzs[warning],  significance(result['Degraded']['zz'], result['Master']['zz'], result['Master']['zerr'])[zzmask][warning], ' x', c='r', markersize=3, label=labels[0])
     pl.plot(truzs[~warning], significance(result['Degraded']['zz'], result['Master']['zz'], result['Master']['zerr'])[zzmask][~warning], 'x', c='k', markersize=3)
     
     pl.plot(_['TRUEZ'][select & np.array(is_lrg)], _['SIG'][select & np.array(is_lrg)], 'x', c='g',    markersize=3, label=labels[1])
     pl.plot(_['TRUEZ'][select & np.array(is_elg)], _['SIG'][select & np.array(is_elg)], 'x', c='b',    markersize=3, label=labels[2])
     pl.plot(_['TRUEZ'][select & np.array(is_qso)], _['SIG'][select & np.array(is_qso)], 'x', c='gold', markersize=3, label=labels[3])

     labels     = [''] * 4
    
  except:
    continue  

  results.append(result)

#
pl.axhline(y=0.0, xmin=0, xmax=1, c='k')
pl.axhline(y=1.0, xmin=0, xmax=1, c='k')
  
pl.xlim(-0.02, 2.1)
pl.ylim(-2.00, 5.0)

pl.xlabel(r'$z_{\rm{True}}$')
pl.ylabel(r"$\log_{10}(|z' - z| \ / \ z_{\rm{err}})$")

pl.legend(loc=1, frameon=True)

ax = pl.gca()

ax.grid(False)

plt.tight_layout()

pl.savefig('zeff.pdf')

print('\n\nDone.\n\n')
