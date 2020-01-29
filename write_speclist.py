import glob


root = '/global/projecta/projectdirs/desi/datachallenge/redwood/spectro/redux/redwood/spectra-64/'

dirs = glob.glob(root + '*')

print(dirs)

ts   = []
es   = []

for _ in dirs:
 files = glob.glob(_ + '/*')    
 ts   += files
 
for _ in ts:
 exp   = _.split('/')[-1]
 es.append(exp)

# write redwood spectra list
with open('redwood_spectra_files.txt', 'w') as filehandle:
    for listitem in ts:
        filehandle.write('%s\n' % listitem)

print('\n\n{} spectra   in redwood simulation'.format(len(ts)))
print('\n\n{} exposures in redwood simulation'.format(len(es)))

print('\n\nDone.\n\n')
