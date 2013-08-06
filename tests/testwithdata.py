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
nf = 1e3
h0 = 2

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

methods = ['GR', 'G4v']

s = srch.inject('GR', [1.0, 3.0, 9.0])

print s

# plt.plot(h, s)
# plt.title('Test')
# plt.savefig('plots/test.png')
# plt.close()
# # 

# a, cov = srch.chi('GR')
# print cov
# exit()
# print time()-start

# 
# 
# hdf = pd.DataFrame(h)
# sdf = pd.DataFrame(s)
# 
# sid, ap = plt.plot(hdf, sdf, '+')
# plt.setp(ap, c='b', label='AP')
# plt.setp(sid, c='c', label='Sid')
# plt.legend()
# plt.legend(numpoints=1)
# 
# plt.savefig('plots/hsAPsidLHO.png')
# plt.close()
# # print cov['f0']
