print 'Loading modules...'

import process
import paths

reload(paths)
reload(process)

paths.makestructure()


# SET PARAMETERS
nf = 2000
ninj = 100

detname = 'H1'
psr = 'J0534+2200'

injection_kinds = ['GR', 'G4v']

search_methods = ['GR', 'G4v']

pd = [0.0, 'p', 'm']

plots = ['hinjrec', 'hinjs', 'hinjlins']


def crab(range=''):

    for kind in injection_kinds:

        for p in pd:
    
            ij = process.InjSearch(detname, psr, nf, kind, p, ninj, rangeparam=[range])
        
            ij.analyze(search_methods)
        
            for plot in plots:
                ij.results.plots(crab, plot, extra_name='S6_range'+range)
            ij.results.save(extra_name=range)
    
