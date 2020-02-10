import os
import glob


# S/N degraded by f * f. 
run   = 'degraded'
files = glob.glob('/global/cscratch1/sd/mjwilson/desi/bluebump/degraded/*')

# Remove redwood output dir. from production list.  
files.remove('/global/cscratch1/sd/mjwilson/desi/bluebump/degraded/redrock')

# Throughput with dip.
# run   = 'thrudip'
# files = glob.glob('/global/cscratch1/sd/mjwilson/desi/bluebump/thrudip/*')
# files.remove('/global/cscratch1/sd/mjwilson/desi/bluebump/thrudip/redrock')

# Redwood rerun with redrock master branch.
# run   = 'redwood'
# file  = open('redwood_spectra_files.txt', mode = 'r')
# files = file.readlines()
# files = [x.replace('\n', '') for x in files]
# file.close()

for i, x in enumerate(files):
  # -- Master redrock run of original redwood files --
  # _ = x.split('/')[-1]
  # x = x + '/spectra-64-{}.fits'.format(_)

  # Degraded or thrudip. 
  _   = x.split('/')[-1].split('-')[2].replace('.fits', '')

  # Only reduce degraded spectra for which the new master exists (i.e. no throughput change) for comparison.
  if os.path.exists(os.environ['CSCRATCH'] + '/desi/bluebump/redwood/redrock/redrock-64-{}.h5'.format(_)):
    ofile = os.environ['CSCRATCH'] + '/desi/bluebump/{}/redrock/redrock-64-{}.h5'.format(run, _)
    zbest = os.environ['CSCRATCH'] + '/desi/bluebump/{}/redrock/zbest-64-{}.fits'.format(run, _)
  
    if os.path.exists(ofile):
      print('{} exists.'.format(ofile))

    else:
      cmd   = 'srun -N 8 -n 256 -c 2 rrdesi_mpi ' + x + ' -o ' + ofile + ' -z ' + zbest

      print(cmd)

      os.system(cmd)
