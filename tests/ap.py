import numpy as np
import pandas as pd
from time import time

from templates import sidereal as sd
from analysis import search
from templates import detector
from templates import source
from templates import antennapatterns as aP
from templates import temps

reload(search)
reload(detector)
reload(source)
reload(aP)

days = 1
t = [x+ 1056563236 for x in range(0, int(days*sd.ss), 60)]

det = 'LHO'
crab = 'J0534+2200'

print "Making up data...",
start=time()
psrlist = source.getpsrlist()
d = np.random.random(len(t)*len(psrlist)).reshape(len(t), len(psrlist))
data = pd.DataFrame(d, index=t, columns=psrlist)
print "Data took %r" % str(start-time())

print "Getting sigma...",
start = time()
sigma = search.segmentsigma(data)
print "Sigma took %r" % str(start-time())

# p = aP.getAP(det, data.index)
# 
# print p['pl']
# 
# 
iota = np.radians(62)

start=time()
a, covmat = search.chiSid(det, iota, data[crab], sigma[crab])
print "Search took %r" % str(start-time())
