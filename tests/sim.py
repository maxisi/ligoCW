import numpy as np
import pandas as pd
from time import time

from templates import sidereal as sd
from analysis import search
from templates import detector
from templates import source
from templates import antennapatterns as aP
from templates import simulate
from analysis import heterodyne, background
import matplotlib.pyplot as plt

reload(search)
reload(detector)
reload(source)
reload(aP)
reload(heterodyne)
reload(background)
reload(simulate)

# PARAMETERS TO SET:
days = 1
nf = 1e5

h0 = 1e-24

####
f = np.linspace(1e-5, 1e-7, nf)

t = [x+ 1056563236 for x in range(0, int(days*sd.ss), 60)]

det = 'LHO'
crab = 'J0534+2200'

psrlist = source.getpsrlist()
d = np.random.random(len(t)*len(psrlist)).reshape(len(t), len(psrlist))
data = pd.DataFrame(d, index=t, columns=psrlist)
####

sign = simulate.Signal('GR', det, h0, t)

a =sign.simulate()[crab]

print type(a)
# print '\n'
# print sign.F['pl']
# print '\n'
# print sign.h