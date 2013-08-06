import numpy as np
import pandas as pd
from time import time

from templates import sidereal as sd
from analysis import search
from templates import detector
from templates import source
from templates import antennapatterns as aP
from templates import temps
from analysis import heterodyne

reload(search)
reload(detector)
reload(source)
reload(aP)
reload(heterodyne)

# PARAMETERS TO SET:
days = 1
nf = 1e5

####
f = np.linspace(1e-5, 1e-7, nf)

t = [x+ 1056563236 for x in range(0, int(days*sd.ss), 60)]

det = 'LHO'
crab = 'J0534+2200'

psrlist = source.getpsrlist()
d = np.random.random(len(t)*len(psrlist)).reshape(len(t), len(psrlist))
data = pd.DataFrame(d, index=t, columns=psrlist)
####

print 'Days: %d' % days
print 'Freqs: %d' % nf


start = time()
heterodyne.createbackground(data, f, det)

# rhB_dir = 'files/background/'
# rh = pd.HDFStore(rhB_dir + det + '0')
# rh2 = pd.HDFStore(rhB_dir + det + '1')
# 
# print rh
# print rh2
print time() - start

# background.rhsave(data, f, det)
# 
# a = heterodyne.het(data[crab], f)
# print a[0]