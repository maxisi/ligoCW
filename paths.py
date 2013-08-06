'''
Contains all paths used by different modules, making it easy to make changes.

File structure:
files/
    analysis/
        results/
    data/
        ephemerides/
    background/
    templates/
        antennapatterns/
    remote/
        source/
        ephemerides/
plots/
        
'''
import os

from numpy import pi

leafList = ('analysis/results/', 'data/ephemerides/', 'background/', 'templates/antennapatterns')
def makestructure():
    print 'Creating file structure:'
    # create files/
    for dir in leafList:
        try:
            os.makedirs('files/' + dir)
            print dir
        except OSError:
            print 'Error: files/%(dir)s already exists' % locals()
    # create plots/
    try:
        os.mkdir('plots/')
        print 'plots/'
    except OSError:
        print 'Error: plots/ already exists'


# Antenna patterns
psrlist = 'files/templates/psrlist.txt'         # list of pulsars
ap = 'files/templates/antennapatterns/'         # antenna patterns

# Data
rhB = 'files/background/'                       # background heterodynes
originalData = 'files/remote/source' # '../../matthew/analyses/S6_all/results'
importedData = 'files/data/dataPitkin_'         # original data in DataFrame

# Detector
vectors = 'files/templates/vectors/'            # detector and wave vectors

# Source
textfromATNF = 'files/templates/psrcat.txt'     # download from ATNF
psrcat = 'files/templates/pulsar_catalogue'     # pulsar data in DF
psrextra = 'files/psrextra.txt'                 # polarization and inclination angles
psrlist = 'files/templates/psrlist.txt'         # names of PSR names to analyze

# Analysis
sigma = 'files/analysis/psrSegmentSigma'        # day-long standard deviation
results = 'files/analysis/results/'             # background h_rec and significance
tmat = 'files/analysis/tmat'                    # basis transformation matrix

# Ephemerides
eph = 'files/remote/ephemerides/' #'usr/share/lalpulsar/'
eph_local = 'files/data/ephemerides/'           # imported ephemerides in HDF5 format


def pname(pinj):
    if pinj== pi/2:
        return 'p'
    elif pinj==-pi/2:
        return 'm'
    elif pinj==0:
        return '0'
    else:
        return str(pinj)