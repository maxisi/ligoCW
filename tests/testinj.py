import numpy as np
import pandas as pd
from time import time
from sys import exit
import matplotlib.pyplot as plt

from templates import sidereal as sd
from analysis import search
from templates import detector
from templates import source
from templates import antennapatterns as aP
from templates import temps
from analysis import data
from analysis import process
import paths
from templates import simulate

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
ninj = 5

hinjrange=[1.e-20, 30]
hinj = np.linspace(hinjrange[0], hinjrange[1], ninj)

####
det = 'LHO'
detector = 'H1'
crab = 'J0534+2200'
other = 'J0023+09'

back_file = pd.HDFStore('files/background/LHO0', 'r')
srch = search.Search(det, crab, back_file[crab])

inj0 = srch.inject('GR', hinj)
plt.plot(srch.data['f0'][:10000])
plt.savefig('test_injout3.png')
print 'plotted' 
plt.close()

# h[psr], s[psr] = injsrch0.generalChi(methods)



