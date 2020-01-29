import numpy as np
import pylab as pl

import astropy.io.fits as fits

from   astropy.table import Table


dat = Table(fits.open('data/throughput/thru-b.fits')[1].data)

