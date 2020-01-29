import os


file  = open('redwood_truth_spectra_files.txt', mode = 'r')
lines = file.readlines()
file.close()

lines = [x.replace('\n','') for x in lines]
exps  = [x.split('/')[-1] for x in lines]
lines = [x + '/truth-64-{}.fits'.format(x.split('/')[-1]) for x in lines]

for i, x in enumerate(lines):
  cmd = 'quickspectra -i ' + x + ' -o ' + os.environ['CSCRATCH'] + '/desi/bluebump/spectra-64-{}.fits'.format(exps[i])

  print(cmd)

  os.system(cmd)
