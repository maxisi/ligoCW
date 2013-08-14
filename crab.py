print 'Loading modules...'

import process
import paths


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


def crab(range='')

    if range=='':
        range = []
        
    for kind in injection_kinds:

        for p in pd:
    
            ij = process.InjSearch(detector, crab, nf, kind, p, ninj, rangeparam=[range])
        
            ij.analyze(search_methods)
        
            for plot in plots:
                ij.results.plots(crab, plot, extra_name='S5_range'+range)
            ij.results.save(extra_name=range)
    
