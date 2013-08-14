#!/usr/bin/env python

print 'Loading modules...'

import process
import paths
import sys


paths.makestructure()


# SET PARAMETERS
nf = 2000
ninj = 100

detector = 'H1'
crab = 'J0534+2200'

injection_kinds = ['GR', 'G4v']

search_methods = ['GR', 'G4v']

pd = [0.0, 'p', 'm']

plots = ['hinjrec', 'hinjs', 'hinjlins']

range = sys.argv
del range[0]

if 'all' in range:
    range = 'all'



# PROCESS
for kind in injection_kinds:

    for p in pd:
    
        ij = process.InjSearch(detector, crab, ninj, kind, pd, ninj)
        
        ij.analyze(search_methods)
        
        ij.results.save(extra_name=range)
    
        for plot in plots:
            ij.results.plots(crab, plot, extra_name='S5_range'+range)