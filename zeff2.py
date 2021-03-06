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

# Read the truth table.
truth   = Table.read(truth, format='fits')

tids    = truth['TARGETID']
tzs     = truth['TRUEZ']

# Remove unwanted columns for e.g. easy viewing. 
del  truth['MOCKID']
del  truth['SEED']
del  truth['TEFF']
del  truth['LOGG']
del  truth['FEH']
del  truth['CONTAM_TARGET']
del  truth['TEMPLATEID']
del  truth['VDISP']

# Nanomaggies to mags.
for x in ['FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2']:
  ff       = truth[x]  
  truth[x] = 22.5 - 2.5 * np.log10(ff)
  
print('\n')

# Add background redshifts affected for each line in output pdf. 
for line, color, label in zip([1549.5, 1908.7, 2799.117, 3727., 3889.0, 4072.3, 4364.436], ['k', 'g', 'chocolate', 'c', 'r', 'm'], ['CIV', 'CIII', 'MgII', 'OII', 'HeI', 'SII', 'OIII']):
  lo = 4300. / line - 1.0
  hi = 4500. / line - 1.0
  
  pl.fill_between(np.arange(lo, hi, 0.01), -2., 5., color=color, alpha=0.2)  

  plt.text((hi + lo) / 2., 4.7, label, {'ha': 'center'})
  
  print('{:.3f}  {:.3f}  {:.3f}'.format(line, lo, hi))

##
tiles      = np.loadtxt('/global/homes/m/mjwilson/desi/bluebump/__files__/redwood_exposure_nums.txt', dtype=str)
labels     = ['BAD Z', 'LRG', 'ELG', 'QSO'] 


results    = []

# Loop over files. 
for i, tile in enumerate(tiles): 
  # Original redwood run files.
  # redwood  = root + '/spectro/redux/redwood/spectra-64/{}/{}/zbest-64-{}.fits'.format(np.int(np.float(tile) / 100.), tile, tile)

  # Rerun of redwood spectra with redrock master: 02/07/20.
  master   = os.environ['CSCRATCH'] + '/desi/bluebump/redwood/redrock/zbest-64-{}.fits'.format(tile)

  # S/N degraded according to throughput via IVAR; redrock master.  
  degrade  = os.environ['CSCRATCH'] + '/desi/bluebump/degraded/redrock/zbest-64-{}.fits'.format(tile) 

  # Run in which actual throughput is degraded for quickspectra;  NOTE:  had no dynamic exposure times so not directly comparable to redwood. 
  # dipthru = os.environ['CSCRATCH'] + '/desi/bluebump/thrudip/redrock/zbest-64-{}.fits'.format(tile)

  # Corresponding Truth file.
  # ttargets = root + '/targets/{}/{}/truth-64-{}.fits'.format(np.int(np.float(tile) / 100.), tile, tile)
  
  result = {} 
  
  # 'Dip', 'Redwood'
  for x, label in zip([master, degrade], ['Master', 'Degraded']):
    # Load a spectroscopic tile.
    result[label] = {}

    zbest         = fits.open(x)[1]
    ztids         = zbest.data['TARGETID']
    
    '''
    isin          = np.ones_like(ztids).astype(bool)

    if x == dipthru:
      ztids       = fits.open(ttargets)[1].data['TARGETID'][:5000]
      
      assert  np.all(ztids == np.sort(ztids))

      matched     = np.sort(result['Redwood']['ztids'][np.isin(result['Redwood']['ztids'], ztids)])
      isin        = np.isin(ztids, matched)

      result[label]['matched'] = np.isin(result['Redwood']['ztids'], matched)      
    '''

    #
    zz            = zbest.data['Z']
    zerr          = zbest.data['ZERR']
    zwarn         = zbest.data['ZWARN']
    
    # TARGETIDS is the spec. tile.
    ztids         = ztids

    # Match TIDS in spec. tile with target ids in truth table with a sorted search.
    index         = np.argsort(tids)
    sorted_tids   = tids[index]
    sorted_index  = np.searchsorted(sorted_tids, ztids)

    yindex        = np.take(index, sorted_index, mode="clip")
    mask          = tids[yindex] != ztids

    trusort       = np.ma.array(yindex, mask=mask)

    # Spec. matched TARGETIDs in truth. 
    _             = tids[trusort][~trusort.mask]

    # Corresponding redshifts.
    truzs         =  tzs[trusort][~trusort.mask]

    # Objects that weren't in the truth tables.  Standards?
    zzmask        = np.isin(ztids, _)

    # for i, x in enumerate(ztids[zzmask]):
    #   print(x, _[i])
    
    # Check that TID matches between truth table and spec. tile was successful.
    assert  np.all(ztids[zzmask] == _)

    # Save to a dict.
    result[label]['zbest']   = zbest
    result[label]['ztids']   = ztids[zzmask]
    result[label]['zz']      = zz[zzmask]
    result[label]['zerr']    = zerr[zzmask]
    result[label]['zwarn']   = zwarn[zzmask] > 0  # Retain only if there was a warning.
    result[label]['zzmask']  = zzmask
    result[label]['tids']    = _
    result[label]['truzs']   = truzs
    result[label]['zbias']   = significance(truzs, zz[zzmask], zerr[zzmask])

    result[label]['nmatch']  = len(ztids)

    ##
    print('\n\nHealpixel:  {}'.format(tile))

    print('Number of targets: {}'.format(len(result[label]['ztids'])))
    print('Number of warnings for {}: {}'.format(label, np.count_nonzero(zwarn)))
    
  # Check that no funny sorting happened between the Master and Degraded runs. 
  assert  np.all(result['Degraded']['ztids'] == result['Master']['ztids'])

  #  Construct the truth table for the targets matched to this spec. tile.
  isin             =  np.isin(truth['TARGETID'], result['Master']['ztids'])

  _                =  Table(truth[isin], copy=True)
  _.sort('TARGETID')

  # Args. to sort the spec. tile by targetid.
  inds             =  np.argsort(result['Master']['ztids'])
  
  # Double check right matching.  
  assert  np.all(_['TARGETID'] == result['Master']['ztids'][inds])
  
  # Calculate the significance of the redshift shift between Degraded and Master.
  _['SIG']         =  significance(result['Degraded']['zz'], result['Master']['zz'], result['Master']['zerr'])[inds]

  # Updates of Master and Degraded results.
  _['MASTERZ']     =  result['Master']['zz'][inds]
  _['DEGRDEZ']     =  result['Degraded']['zz'][inds]

  _['MASTERZERR']  =  result['Master']['zerr'][inds]
  _['DEGRDEZERR']  =  result['Degraded']['zerr'][inds]

  _['MASTERZWARN'] =  result['Master']['zwarn'][inds]
  _['DEGRDEZWARN'] =  result['Degraded']['zwarn'][inds]
  
  # Remove white space.                                                                                                                                                                                              
  _['TEMPLATETYPE'] = [x.strip() for x in _['TEMPLATETYPE']]
  
  # Define the reference sample. 
  _['IN_MASTER']       = (_['MASTERZWARN'] == 0) & (_['TRUEZ'] <= 2.1)
  _['IN_DEGRDE']       = (_['DEGRDEZWARN'] == 0) & (_['TRUEZ'] <= 2.1)
  
  # Remember: SIG is log10(significance)!
  _['NATURALV_MASTER'] = np.abs(_['MASTERZ'] - _['TRUEZ']) / (1. + _['TRUEZ'])
  _['NATURALV_DEGRDE'] = np.abs(_['DEGRDEZ'] - _['TRUEZ']) / (1. + _['TRUEZ'])

  ##  1000. km/s.
  _['CATFAIL_MASTER']  = _['IN_MASTER'] & (_['NATURALV_MASTER'] > 0.003)
  _['CATFAIL_DEGRDE']  = _['IN_DEGRDE'] & (_['NATURALV_DEGRDE'] > 0.003)
  
  print()
  print('Solved for: {}'.format(i))
  print(_)
  
  is_elg = [x == 'ELG' for x in _['TEMPLATETYPE']]
  is_lrg = [x == 'LRG' for x in _['TEMPLATETYPE']]
  is_qso = [x == 'QSO' for x in _['TEMPLATETYPE']]

  pl.plot(_['TRUEZ'][(_['MASTERZWARN'] == 0) & is_lrg], (_['DEGRDEZERR'] / _['MASTERZERR'])[(_['MASTERZWARN'] == 0) & is_lrg], marker='x', c='r',    markersize=3, label='', lw=0)
  pl.plot(_['TRUEZ'][(_['MASTERZWARN'] == 0) & is_elg], (_['DEGRDEZERR'] / _['MASTERZERR'])[(_['MASTERZWARN'] == 0) & is_elg], marker='x', c='b',    markersize=3, label='', lw=0)  
  pl.plot(_['TRUEZ'][(_['MASTERZWARN'] == 0) & is_qso], (_['DEGRDEZERR'] / _['MASTERZERR'])[(_['MASTERZWARN'] == 0) & is_qso], marker='x', c='gold', markersize=3, label='', lw=0)

  results.append(_)

  #if i > 20:
  #  break
  
##  Sample stats.
##  -------------------------
##  Counts of the simulated galaxies in each target class for MASTERZWARN == 0.   
m_nbgs    = 0
m_nlrg    = 0
m_nelg    = 0
m_nqso    = 0
m_nstar   = 0

##  Counts of the simulated galaxies in each target class for DEGRADEDZWARN == 0.                                                                                                                                                 
d_nbgs    = 0
d_nlrg    = 0
d_nelg    = 0
d_nqso    = 0
d_nstar   = 0

##  Catastrophic failures with respect to the MASTER z. 
cm_nbgs   = 0
cm_nlrg   = 0
cm_nelg   = 0
cm_nqso   = 0
cm_nstar  = 0

##  Catastrophic failures with respect to the DEGRADED z.                                                                                                                                                                        
cd_nbgs   = 0
cd_nlrg   = 0
cd_nelg   = 0
cd_nqso   = 0
cd_nstar  = 0

##  Check on all available TEMPLATETYPES.
##
##  TEMPLATETYPE
##  ------------
##       ELG
##       LRG
##       QSO
##       STAR

utemps, cnts = np.unique(_['TEMPLATETYPE'], return_counts=True)

print('\n\n')
print(utemps)
print(cnts)
print('\n\n')

for i, _ in enumerate(results):
  m_nbgs   += np.count_nonzero((_['IN_MASTER']      & (_['TEMPLATETYPE'] == 'BGS')))
  m_nlrg   += np.count_nonzero((_['IN_MASTER']      & (_['TEMPLATETYPE'] == 'LRG')))
  m_nelg   += np.count_nonzero((_['IN_MASTER']      & (_['TEMPLATETYPE'] == 'ELG')))
  m_nqso   += np.count_nonzero((_['IN_MASTER']      & (_['TEMPLATETYPE'] == 'QSO')))
  m_nstar  += np.count_nonzero((_['IN_MASTER']      & (_['TEMPLATETYPE'] == 'STAR')))

  d_nbgs   += np.count_nonzero((_['IN_DEGRDE']      & (_['TEMPLATETYPE'] == 'BGS')))
  d_nlrg   += np.count_nonzero((_['IN_DEGRDE']      & (_['TEMPLATETYPE'] == 'LRG')))
  d_nelg   += np.count_nonzero((_['IN_DEGRDE']      & (_['TEMPLATETYPE'] == 'ELG')))
  d_nqso   += np.count_nonzero((_['IN_DEGRDE']      & (_['TEMPLATETYPE'] == 'QSO')))
  d_nstar  += np.count_nonzero((_['IN_DEGRDE']      & (_['TEMPLATETYPE'] == 'STAR')))
  
  cm_nbgs  += np.count_nonzero((_['CATFAIL_MASTER'] & (_['TEMPLATETYPE'] == 'BGS')))
  cm_nlrg  += np.count_nonzero((_['CATFAIL_MASTER'] & (_['TEMPLATETYPE'] == 'LRG')))
  cm_nelg  += np.count_nonzero((_['CATFAIL_MASTER'] & (_['TEMPLATETYPE'] == 'ELG')))
  cm_nqso  += np.count_nonzero((_['CATFAIL_MASTER'] & (_['TEMPLATETYPE'] == 'QSO')))
  cm_nstar += np.count_nonzero((_['CATFAIL_MASTER'] & (_['TEMPLATETYPE'] == 'STAR')))

  cd_nbgs  += np.count_nonzero((_['CATFAIL_DEGRDE'] & (_['TEMPLATETYPE'] == 'BGS')))
  cd_nlrg  += np.count_nonzero((_['CATFAIL_DEGRDE'] & (_['TEMPLATETYPE'] == 'LRG')))
  cd_nelg  += np.count_nonzero((_['CATFAIL_DEGRDE'] & (_['TEMPLATETYPE'] == 'ELG')))
  cd_nqso  += np.count_nonzero((_['CATFAIL_DEGRDE'] & (_['TEMPLATETYPE'] == 'QSO')))
  cd_nstar += np.count_nonzero((_['CATFAIL_DEGRDE'] & (_['TEMPLATETYPE'] == 'STAR')))
  
## 
print('\n\nSummary stats:')
print('# galaxies in z <= 2.1 sample (no MASTER ZWARN).')
print()
print('# BGSs  in sample: {}'.format(m_nbgs))
print('# LRGs  in sample: {}'.format(m_nlrg))
print('# ELGs  in sample: {}'.format(m_nelg))
print('# QSOs  in sample: {}'.format(m_nqso))
print('# STARs in sample: {}'.format(m_nstar))
print()
print('# galaxies in z <= 2.1 sample (no DEGRADED ZWARN).')
print('# BGSs  in sample: {}'.format(d_nbgs))
print('# LRGs  in sample: {}'.format(d_nlrg))
print('# ELGs  in sample: {}'.format(d_nelg))
print('# QSOs  in sample: {}'.format(d_nqso))
print('# STARs in sample: {}'.format(d_nstar))
print()
print('# BGSs  master cat. fails: {}'.format(cm_nbgs))
print('# LRGs  master cat. fails: {}'.format(cm_nlrg))      
print('# ELGs  master cat. fails: {}'.format(cm_nelg))
print('# QSOs  master cat. fails: {}'.format(cm_nqso))
print('# STARs master cat. fails: {}'.format(cm_nstar))
print()
print('# BGSs  degraded cat. fails: {}'.format(cd_nbgs))
print('# LRGs  degraded cat. fails: {}'.format(cd_nlrg))
print('# ELGs  degraded cat. fails: {}'.format(cd_nelg))
print('# QSOs  degraded cat. fails: {}'.format(cd_nqso))
print('# STARs degraded cat. fails: {}'.format(cd_nstar))

##
pl.axhline(y=0.0, xmin=0, xmax=1, c='k')
  
pl.xlim(-0.02, 2.1)
pl.ylim( 0.00, 3.0)

pl.xlabel(r'$z_{\rm{True}}$')

# pl.ylabel(r"$\log_{10}(|z' - z| \ / \ z_{\rm{err}})$")
pl.ylabel(r"$(\sigma{_z}' / \sigma{_z})$")


pl.legend(loc=1, frameon=True)

ax = pl.gca()

ax.grid(False)

plt.tight_layout()

pl.savefig('plots/zeff2.pdf')

##
pl.clf()

print('\n\nDone.\n\n')
