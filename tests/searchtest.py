import numpy as np
import pandas as pd
from time import time

from templates import sidereal as sd
from analysis import search
from templates import detector
from templates import source
from templates import antennapatterns as aP
from templates import temps
from analysis import heterodyne, background
from templates import simulate

reload(search)
reload(detector)
reload(source)
reload(aP)
reload(heterodyne)
reload(background)
reload(temps)
reload(simulate)

# PARAMETERS TO SET:
days = 1
nf = 1e5

####
f = np.linspace(1e-5, 1e-7, nf)

t = [x+ 1056563236 for x in range(0, int(days*sd.ss), 60)]

det = 'LHO'
crab = 'J0534+2200'
incl = np.radians(62)

psrlist = source.getpsrlist()
d = np.random.random(len(t)*len(psrlist)).reshape(len(t), len(psrlist))
data = pd.DataFrame(d, index=t, columns=psrlist)
sigma = search.segmentsigma(data)
####

a, cov = search.chi('GR', det, data, sigma)
a2, cov2 = search.chi('G4v', det, data, sigma)

print a[crab]
print a2 [crab]
print '\n'
print cov[crab]
print cov2[crab]