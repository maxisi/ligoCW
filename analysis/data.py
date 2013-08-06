import numpy as np
import pandas as pd
import tables as tb
from sys import exit

from templates.source import getpsrlist
import paths

threshold = 100 # max in hdf5: 1e4
print 'Warning: adjust threshold. Data line 9'

def het(vector, f, *arg):
    '''
    Heterodynes vector at frequencies f. Preferred input for vector is series indexed
    over time; f can be a list or an array. Returns DataFrame.
    '''
    print 'Ready to heterodyne.',
    
    if len(arg)==0:
        try:
            t = vector.index.tolist()
        except AttributeError:
            print 'ERROR: no time vector for heterodyne.'
            exit(0)
            
    elif len(arg)==1:
        t = arg[0]
        
    else:
        print 'ERROR: het needs input time or indexed vector, not %d extra arguments.' % len(arg)
        exit(0)
    
    temp = np.exp(2*np.pi*1j*np.multiply.outer(f, t))
    print 'Template created.'
    
    try:
        template = pd.DataFrame(temp, index=['f' + str(x) for x in range(0,len(f))], columns=t)
    except ValueError:
        template = pd.Series(temp, index=t)
    
    rh = vector*template
    return rh.T
    

def background(data, freq, det):
    '''
    Re heterodynes and saves data at frequencies f. Number of heterodynes is determined by
    f and data can be for more than one pulsar
    '''
    nsets = int(len(freq)/threshold)
    print '\tWill create background in %(nsets)d files of 1e4 instantiations each.' % locals()

    if nsets<1:
        nsets = 1
    else:
        pass
    
    fset = {n : freq[n*threshold:min(len(freq),(n+1)*threshold)] for n in range(0, nsets)}
    path = [paths.rhB + 'back_' +  det + str(n) for n in range(0, nsets)]

    for n in range(0, nsets):
        rh = pd.HDFStore(path[n])
        
        try:
            for psr in [s.strip('/') for s in data.keys()]:
                print psr
                rh[psr] = het(data[psr], fset[n])
                
        except AttributeError:
            print "Atribute Error in data.background"
            rh[data.name] = het(data, fset[n])
        rh.close()
    

def imp(detector):
    '''
    Return DF with original data (col: PSR; index: t). Assuming execution on ATLAS.
    '''
        
    psrlist = getpsrlist()
    try:
        d = pd.HDFStore(paths.importedData + detector + '.hdf5')

        for psr in psrlist:
            psrfile = 'finehet_' + psr + '_' + detector
        
            try:
                psrpath = paths.originalData + '/data' + detector + '/' + psrfile
                dParts = pd.read_table(psrpath, sep='\s+', names= [None, 'Re', 'Im'], header=None, index_col=0)
                
            except IOError:
            
                try:
                    psrpath = paths.originalData + '/' + psr + '_' + detector + '/data' + detector + '/' + psrfile
                    dParts = pd.read_table(psrpath, sep='\s+', names= [None, 'Re', 'Im'], header=None, index_col=0)
                    
                except IOError:
                    print 'Cannot find %s data for PSR %s in %s.' % (detector, psr, psrpath),
                    print 'Should I...'
                    print '\t1. Provide path\n\t2. Abort'
                    opt = raw_input('>')
                    
                    if opt==1:
                    
                        try:
                            psrpath = raw_input('Enter path: ')
                            dParts = pd.read_table(psrpath, sep='\s+', names= [None, 'Re', 'Im'], header=None, index_col=0)
                        
                        except IOError:
                            print "File not found. One more try?"
                            psrpath = raw_input('Enter path: ')
                            dParts = pd.read_table(psrpath, sep='\s+', names= [None, 'Re', 'Im'], header=None, index_col=0)
                    
                    else:
                        print 'Could not find %(detector)s data for PSR %(psr)s.' % locals()
                        print 'Exiting at analysis/data ln 108'
                        exit()

            dWhole = dParts['Re']+dParts['Im']*1j
            d[psr] = dWhole
            
    finally:
        d.close()
    
    
def get(detector, psr=None):
    '''
    Retrieves original heterodyned data for pulsars in list.
    Imports data from M.Pitkin if necessary.
    '''
        
    print 'Getting %(detector)s original heterodyned %(detector)s data.' % locals()

    try:
        d = pd.HDFStore(paths.importedData + detector + '.hdf5')
        psrlist = getpsrlist()
        
        if set(psrlist)!=set([s.strip('/') for s in d.keys()]):
            print 'Not all pulsars are present in file.',
#             print 'PSR list: %r' % psrlist
#             print 'File list: %r' % [s.strip('/') for s in d.keys()]
            print 'Refreshing...',
            d.close()
            imp(detector)
            d = pd.HDFStore(paths.importedData + detector + '.hdf5')
            
            if psr == None:
                print 'Success.'
                return d
            else:
                dpsr = d[psr]
                d.close()
                print 'Success.'
                return dpsr
        else:
            print 'File is ready to go.',
            
            if psr == None:
                print 'Success.'
                return d
            else:
                dpsr = d[psr]
                d.close()
                print 'Success.'
                return dpsr
                
    except KeyError:
        print 'No record for %(detector)s found. Importing it from M.Pitkin...' % locals(),
    
        imp(detector)
        d = pd.HDFStore(paths.importedData + detector + '.hdf5')
        if psr == None:
            print 'Success.'
            return d
        else:
            dpsr = d[psr]
            d.close()
            print 'Success.'
            return dpsr