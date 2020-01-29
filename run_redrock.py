import os
import glob


# S/N degraded. 
run   = 'degraded'
files = glob.glob('/global/cscratch1/sd/mjwilson/desi/bluebump/degraded/*')
files.remove('/global/cscratch1/sd/mjwilson/desi/bluebump/degraded/redrock')

# Redwood rerun.
# run   = 'redwood'
# file  = open('redwood_spectra_files.txt', mode = 'r')
# files = file.readlines()
# files = [x.replace('\n', '') for x in files]
# file.close()

for i, x in enumerate(files):
  #_     = x.split('/')[-1]
  #x     = x + '/spectra-64-{}.fits'.format(_)

  _     = x.split('/')[-1].split('-')[2].replace('.fits', '')

  ofile = os.environ['CSCRATCH'] + '/desi/bluebump/{}/redrock/redrock-64-{}.h5'.format(run, _)
  zbest = os.environ['CSCRATCH'] + '/desi/bluebump/{}/redrock/zbest-64-{}.fits'.format(run, _)
  
  if os.path.exists(ofile):
    print('{} exists.'.format(ofile))

  else:
    cmd   = 'rrdesi ' + x + ' -o ' + ofile + ' -z ' + zbest

    print(cmd)

    os.system(cmd)
