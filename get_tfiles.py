import glob


root = '/global/projecta/projectdirs/desi/datachallenge/redwood/targets/'

dirs = glob.glob(root + '*')
dirs.remove('/global/projecta/projectdirs/desi/datachallenge/redwood/targets/qa')

print(dirs)

ts   = []
es   = []

for _ in dirs:
 files = glob.glob(_ + '/*')    
 ts   += files
 
for _ in ts:
 exp   = _.split('/')[-1]
 es.append(exp)

# write redwood exposure list
with open('redwood_exposure_nums.txt', 'w') as filehandle:
    for listitem in es:
        filehandle.write('%s\n' % listitem)

# write redwood truth spectra list
with open('redwood_truth_spectra_files.txt', 'w') as filehandle:
    for listitem in ts:
        filehandle.write('%s\n' % listitem)

print('\n\n{} truth spectra in redwood simulation'.format(len(ts)))
print('\n\n{} exposures     in redwood simulation'.format(len(es)))

print('\n\nDone.\n\n')
