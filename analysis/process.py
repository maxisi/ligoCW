from sys import exit
import matplotlib.pyplot as plt
import os

import numpy as np
import pandas as pd

# from analysis import search
# from templates import simulate
from templates2 import detnames
import paths


class Data(object):

    def __init__(self, detector, psr, finj):
        self.detector = detector
        self.det = detnames(detector)
        self.psr = psr
        self.finj = finj
        
        # background info
        self.backdir = paths.rhB + self.det + '/' + self.psr + '/'
        self.backname = 'back_' + self.psr + '_' + self.det
        
        # data info
        self.datadir = paths.importedData + self.psr + '_' + self.detector + '.hdf5'
        self.finehet = 'finehet_' + self.psr + '_' + self.detector


    def imp(self):
        '''
        Return DF with original data (col: PSR; index: t). Assuming execution on ATLAS.
        '''
        
        struct = '/data' + self.detector + '/' + self.finehet
        pathOptions = [
                     paths.originalData + struct,
                     paths.originalData + '/' + self.psr + '_' + self.detector + struct
                     ]
        
        try:
            d = pd.HDFStore(self.datadir, 'w')
                
            for p in pathOptions:
                try:
                    dParts = pd.read_table(p, sep='\s+', names= [None, 'Re', 'Im'], header=None, index_col=0)
                except IOError:
                    pass
            
            # check file was found
            try:
                dParts
            except NameError:
                print 'Could not find %s data for PSR %s in:\n' % (self.detector, self.psr),
                for p in pathOptions:
                    print '\t%(p)s' % locals()
                    
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
                    print 'Exiting at analysis/process ln 77'
                    exit()

            dWhole = dParts['Re']+dParts['Im']*1j
            d[self.psr] = dWhole
            
        finally:
            d.close()


    def get(self):
        '''
        Retrieves original heterodyned data for pulsars in list.
        Imports data from M.Pitkin if necessary.
        '''
        
        print 'Getting %(detector)s original heterodyned %(detector)s data.' % locals()

        try:
            d = pd.HDFStore(self.datadir)
        
            if self.psrlist not in [s.strip('/') for s in d.keys()]:
                print 'PSR not present in file.',
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
                
                
    def background(self, data):
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
        path = [self.backpath + str(n) for n in range(0, nsets)]

        for n in range(0, nsets):
            
            rh = pd.HDFStore(path[n], 'w')
        
            try:
                for psr in [s.strip('/') for s in data.keys()]:
                    print psr
                    rh[psr] = het(data[psr], fset[n])
                
            except AttributeError:
                print "Atribute Error in data.background"
                rh[data.name] = het(data, fset[n])
            rh.close()    


class InjSearch(object):
    
    def __init__(self, detector, psr, nfreq, injkind, ninj, frange=[1.0e-7, 1.0e-5], hinjrange=[1.0E-27, 1.0E-24]):
        self.detector = detector
        self.det = detnames(detector)
        self.psr = psr
        
        # rehet info
        self.freq = np.linspace(frange[0], frange[1], nfreq)
        
        # injection info
        self.hinj = np.linspace(hinjrange[0], hinjrange[1], ninj)
        self.injkind = injkind  
        
        # background path info
        self.backdir = path.rhB + self.det + '/' + self.psr + '/'

        
        print hinjrange
        
    
    
    def checkfiles(self):
        '''
        Checks background required for search exits and creates it if needed.
        Returns filename list.
        '''
        print 'Checking background files exist...',
        
        files = [name for name in os.listdir(paths.rhB) if self.det in name]
        nfiles = len(files)
        ninst = nfiles * data.threshold # this might not always be true!
    
        if ninst < len(self.freq):
            print 'not enough files for %s. Generating:' %self.det,
            print '(1) Get seed data;',
            d = data.get(self.detector)
            print '(2) Heterodyne.'
            data.background(d, self.freq, self.det)
            print 'Success.'
            d.close()
            
        else:
            print 'they do.\n'
            
        files = [name for name in os.listdir(paths.rhB) if self.det in name]
        
        return files
        
        
    def analyze(self, methods, injpdif=np.pi/2.):
        '''
        Loads a background file, analyzes it feeding it to search by PSR, saves results 
        by pulsar. Loops over files and returns dictionary with h_rec and sig dataframes
        for each PSR (i.e. there are 2 DFs for each PSR key).
        '''
        files = self.checkfiles()
        nfiles = len(files)
        back_paths = [paths.rhB + filename for filename in files]
        print 'Analyzing %d files.' % nfiles
        
        ## Set up injections
        nsets = int(len(self.freq)/data.threshold)
        if nsets<1:
            nsets = 1
        else:
            pass
        hperfile = int(len(self.hinj)/nsets)
        injset = {back_paths[n] : self.hinj[hperfile*n:(n+1)*hperfile] for n in range(0, nsets)}

        ## Process first file to get PSR list and initialize storage variables.
        try:
            file0 = pd.HDFStore(back_paths[0], 'r')
            psrlist = [s.strip('/') for s in file0.keys()]
        
            h = {}
            s = {}
            for psr in psrlist:  
                          
                injsrch0 = search.Search(self.det, psr, file0[psr])
                inj0 = injsrch0.inject(self.injkind, injset[back_paths[0]], pdif=injpdif)
                h[psr], s[psr] = injsrch0.generalChi(methods)
                                
                if psr==psrlist[0]:
                    inj = inj0

        finally:
            file0.close()
            
        ## Loop over rest of the files
        for path in back_paths[1:]:
            try:
                back_file = pd.HDFStore(path, 'r')
            
                for psr in psrlist:
                    injsrch = search.Search(self.det, psr, back_file[psr])
                    inj_file = injsrch.inject(self.injkind, injset[path], pdif=injpdif)
                    
                    h_file, s_file = injsrch.generalChi(methods)

                    h[psr] = pd.concat([h[psr], h_file], ignore_index=True)
                    s[psr] = pd.concat([s[psr], s_file], ignore_index=True)
                    
                    if psr==psrlist[0]:
                        inj = np.concatenate([inj, inj_file])
            finally:
                back_file.close()
                
        ## Save
        hs_path = paths.results + self.detector + '/' + self.detector + '_' + self.injkind + '_' + paths.pname(injpdif) + '_'
        try:
            print 'saving h'
            h_save = pd.HDFStore(hs_path + 'h', 'w')
            for key, psr_h in h.iteritems():
                h_save[key] = psr_h
        except IOError:
            print 'ERROR: %(hs_path)s does not exist.' % locals()
        finally:
            h_save.close()
            
        try:
            s_save = pd.HDFStore(hs_path + 's', 'w')
            for key, psr_s in s.iteritems():
                s_save[key] = psr_s
        finally:
            s_save.close()

        try:
            os.remove(hs_path + 'hinj.npy')
        except OSError:
            pass
            
        np.save(hs_path + 'hinj', inj)

        return h, s, inj
