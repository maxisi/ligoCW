#!/usr/bin/ python

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
nf = 1e4
ninj = 200

####
det = 'LHO'
detector = 'H1'
crab = 'J0534+2200'
other = 'J0023+09'

methods = ['GR', 'G4v']

injsrch = process.InjSearch(detector, nf, 'GR', ninj)

_, _, _ = injsrch.analyze(methods, injpdif=np.pi/2.)


