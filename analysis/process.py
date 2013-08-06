from sys import exit
import matplotlib.pyplot as plt
import os
import datetime

import numpy as np
import pandas as pd

# from analysis import search
# from templates import simulate
from templates2 import detnames, System
import paths


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


class Data(object):

    def __init__(self, detector, psr):
        self.detector = detector
        self.det = detnames(detector)
        self.psr = psr
        
        # data info
        self.datadir = paths.importedData + self.psr + '_' + self.detector + '.hdf5'
        self.seedname = 'finehet_' + self.psr + '_' + self.detector


    def imp(self):
        '''
        Return DF with original data (col: PSR; index: t). Assuming execution on ATLAS.
        '''
        
        struct = '/data' + self.detector + '/' + self.seedname
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

            self.finehet = dParts['Re']+dParts['Im']*1j
            d[self.psr] = self.finehet
            
        finally:
            d.close()


    def get(self):
        '''
        Retrieves original heterodyned data for pulsars in list.
        Imports data from M.Pitkin if necessary.
        '''
        
        try:
            d = pd.HDFStore(self.datadir, 'r')
        
            try:
                self.finehet = d[self.psr]
            except KeyError:
                # file is empty or is corrupted
                self.imp(detector)
        
        except IOError:
            self.imp()
        
        finally:
            d.close()
                

class Background(object):
    def __init__(self, detector, psr, freq, filesize=100):
        # data
        self.seed = Data(detector, psr)
        self.seed.get()
        
        # background info
        self.freq = freq
        self.filesize = filesize      # number of series per file. Adjust!
        
        self.nsets = int(len(freq)/filesize)
        if self.nsets<1:
            self.nsets = 1
            
        self.fset = {n : freq[n*filesize:min(len(freq),(n+1)*filesize)] for n in range(self.nsets)}
        
        # storing info
        self.backdir = paths.rhB + self.seed.det + '/' + psr + '/'
        self.backname = 'back_' + psr + '_' + self.seed.det + '_'
    
             
    def create(self):
        '''
        Re heterodynes and saves data at frequencies f. Number of heterodynes is determined by
        f and data can be for more than one pulsar
        '''
        
        # create background directory
        try:
            os.makedirs(self.backdir)
        except OSError:
            pass
            
        # create background files
        for n in range(self.nsets):
            path = self.backdir + self.backname + str(n)
            
            try:
                rh = pd.HDFStore(path, 'w')
                rh[self.seed.psr] = het(self.seed.finehet, self.fset[n])
            finally:
                rh.close()
        
        # create log
        now = datetime.datetime.now()
        comments = '# ' + self.seed.detector + '\n# ' + self.seed.psr + '\n# ' + str(now) + '\n'
        fileinfo = 'nsets\tfilesize\n' + str(self.nsets) + '\t' + str(self.filesize)
        try:
            f = open(self.backdir + 'log.txt', 'w')
            f.write(comments + fileinfo)
        finally:
            f.close()

    def get(self):
        '''
        Checks background required for search exits and creates it if needed.
        Returns filename list.
        '''
        # read log
        try:
            readme = pd.read_table(self.backdir + 'log.txt', sep='\s+', skiprows=3)
            log_nfiles = readme['nsets'].ix[0]
            log_filesize = readme['filesize'].ix[0]
            log_nfreq = log_nfiles * log_filesize
            
            # get actual number of background files in directory
            files = [name for name in os.listdir(self.backdir) if 'back' in name]
            nfiles = len(files)
            
            if nfiles!=log_nfiles or log_nfreq!=len(self.freq):
                self.create()
        except IOError:
            # no log found
            self.create()
# 
#         files = [name for name in os.listdir(paths.rhB) if self.det in name]
#         
#         return files


class InjSearch(object):
    
    def __init__(self, detector, psr, nfreq, injkind, ninj, frange=[1.0e-7, 1.0e-5], hinjrange=[1.0E-27, 1.0E-24], filesize=100):
        
        self.info = System(detector, psr)
        
        # rehet info
        self.freq = np.linspace(frange[0], frange[1], nfreq)
        
        self.background = Background(detector, psr, self.freq, filesize)
        
        # injection info
        self.hinj = np.linspace(hinjrange[0], hinjrange[1], ninj)
        self.injkind = injkind  
        
        print hinjrange
        
    def analyze(self, methods, injpdif=np.pi/2.):
        '''
        Loads a background file, analyzes it feeding it to search by PSR, saves results 
        by pulsar. Loops over files and returns dictionary with h_rec and sig dataframes
        for each PSR (i.e. there are 2 DFs for each PSR key).
        '''
        
        # background
        self.background.get()
        
        print 'Analyzing %d files.' % self.background.nsets
        
        ## Set up injections
        self.hperfile = int(len(self.hinj)/self.background.nsets)
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
