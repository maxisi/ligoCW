import numpy as np
import pandas as pd
from time import time

from templates import sidereal as sd
from analysis import search
from templates import detector
from templates import source
from templates import antennapatterns as aP
from templates import temps
from templates import simulate
from analysis import data
from analysis import process
import paths

reload(search)
reload(process)
reload(detector)
reload(source)
reload(aP)
reload(data)
reload(paths)
reload(simulate)

# PARAMETERS TO SET:
nf = 100
h0 = [2.0, 3.0]

days = 1
t = [x+ 1056563236 for x in range(0, int(days*sd.ss), 60)]

####
det = 'LHO'
detector = 'H1'
crab = 'J0534+2200'
other = 'J0023+09'

# srch = process.InjSearch(detector, nf, h0)
# 
# n = srch.analyze()
# 
# print n

import matplotlib.pyplot as plt

back_file = pd.HDFStore('files/background/LHO0', 'r')
srch = search.Search(det, crab, back_file[crab])

# a, cov=srch.chi('GR')

sig = simulate.Signal(det, h0, t, crab)

start = time()
s = sig.sim('GR', h0)
print time()-start

print '\n'
print s
