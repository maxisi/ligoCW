from sys import exit
import matplotlib.pyplot as plt
import os
import datetime
from time import time

import numpy as np
import pandas as pd

# from analysis import search
# from templates import simulate
import templates
from templates import detnames, System
from templates_OLD import sidereal as sd
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
        self.dir = paths.rhB + self.seed.det + '/' + psr + '/'
        self.name = 'back_' + psr + '_' + self.seed.det + '_'
        self.path = self.dir + self.name
    
             
    def create(self):
        '''
        Re heterodynes and saves data at frequencies f. Number of heterodynes is determined by
        f and data can be for more than one pulsar
        '''
        
        # create background directory
        try:
            os.makedirs(self.dir)
        except OSError:
            pass
            
        # create background files
        for n in range(self.nsets):
            path = self.dir + self.name + str(n)
            
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
            f = open(self.dir + 'log.txt', 'w')
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
            readme = pd.read_table(self.dir + 'log.txt', sep='\s+', skiprows=3)
            log_nfiles = readme['nsets'].ix[0]
            log_filesize = readme['filesize'].ix[0]
            log_nfreq = log_nfiles * log_filesize
            
            # get actual number of background files in directory
            files = [name for name in os.listdir(self.dir) if 'back' in name]
            nfiles = len(files)
            
            if nfiles!=log_nfiles or log_nfreq!=len(self.freq) or log_filesize!=self.filesize:
                self.create()
        except IOError:
            # no log found
            self.create()
# 
#         files = [name for name in os.listdir(paths.rhB) if self.det in name]
#         
#         return files


class Results(object):
    
    def __init__(self, info, methods, hinj, injkind=None, injpdif=None):
        self.syst = info
        self.methods = methods
        self.hinj = hinj
        self.injkind = injkind
        self.injpdif = injpdif
        
        self.h = pd.DataFrame(columns = methods)
        self.s = pd.DataFrame(columns = methods)
        
        self.dir = paths.results + self.syst.detector + '/' + self.syst.psr + '/' 
        self.name = self.syst.psr + '_' + self.syst.detector + '_' + self.injkind + '_' + paths.pname(injpdif)
        self.path = self.dir + self.name
        self.issaved =  False
        
    def save(self):
    
        try:
            os.makedirs(self.dir)
        except OSError:
            pass
            
        try:
            f = pd.HDFStore(self.path, 'w')
            f['h'] = self.h
            f['s'] = self.s
        finally:
            f.close()
            
        self.issaved = True
        
                    
    def checkInjections(self):
        if self.hinj == h.index:
            print 'All good'
        else:
            print 'Something is wrong...'


class InjSearch(object):
    
    def __init__(self, detector, psr, nfreq, injkind, ninj, frange=[1.0e-7, 1.0e-5], hinjrange=[1.0E-27, 1.0E-24], filesize=100):
        
        self.syst = System(detector, psr)
        
        # rehet info
        self.freq = np.linspace(frange[0], frange[1], nfreq)
        
        self.background = Background(detector, psr, self.freq, filesize)
        
        # injection info
        self.hinj = np.linspace(hinjrange[0], hinjrange[1], ninj)
        self.injkind = injkind
        
        self.hperfile = int(len(self.hinj)/self.background.nsets)
        self.injset = {n : self.hinj[self.hperfile*n:(n+1)*self.hperfile] for n in range(self.background.nsets)}
        
        print hinjrange
        
    def analyze(self, methods, injpdif=np.pi/2.):
        
        # background
        self.background.get()
        print 'Analyzing %d files.' % self.background.nsets
        
        # results
        self.results = Results(self.syst, methods, self.hinj, injkind=self.injkind, injpdif=injpdif)
        
        srch = Search(self.syst, methods=methods)
        
        # loop over files
        for n in range(self.background.nsets):
            try:
                back_file = pd.HDFStore(self.background.path + str(n), 'r')
                
                # initialize data
                srch.initialize_data(back_file[self.syst.psr])
                # inject
                srch.inject(self.syst, self.injkind, self.injset[n], pdif=injpdif)
                
                h_file, s_file = srch.generalChi()
            
                self.results.h.append(h_file)
                self.results.s.append(s_file)
                
            finally:
                back_file.close()
                
        ## Save
        self.results.save()


## SEARCH
def segmentsigma(data):
    '''
    Splits data into day-long segments and returns their standard deviation.
    '''
    # Check orientation    
    t = data.index
    interval_length= sd.ss # ADJUST!
    print 'Taking std over %f second-long intervals:' % interval_length,

    # Slice up data into day-long bins and get groupby stats (see Ch 9 of Python for Data Analysis).
    bins = np.arange(t[0]-interval_length, t[-1]+interval_length, interval_length)
    slices = pd.cut(t, bins, right=False)
    print 'sliced,',

    def getsigma(group):
#         s = np.std(group)
        g = np.array(group.tolist())
        s = np.std(g)
        return s
        #return group.std(ddof=0) # this is pd unbiased 1/(n-1), should use np.std 1/n?
        
    print 'std taken,',
    grouped = data.groupby(slices) # groups by bin
    sigmagroups= grouped.apply(getsigma) # gets std for each bin
    print 'grouped,',

    # Create standard deviation time series 
    s = [sigmagroups.ix[slices.labels[t_index]] for t_index in range(0,len(t)) ]
    sigma = pd.Series(s, index=t)
    print 'done.'
    
    return sigma
    
    
def getsigma(d, psr):
    print 'Retrieving segment standard deviation for PSR %(psr)s...' % locals(),
    try:
        s = pd.HDFStore(paths.sigma)
        try:
            sigma = s[psr]

        except KeyError:
            print 'PSR not in file.'
            sigma = segmentsigma(d)
            s[psr] = sigma
    finally:
        s.close()
        
    print 'Sigma is ready.'
    return sigma
    

def herm(A):
    '''
    Returns hermitian conjugate of A.
    '''
    return np.array(A).conjugate().T
    
    
def svd_analysis(A, b):
    '''
    Performs SVD analysis on design matrix A and data y. "components" is a list with the
    name of the regression components.
    '''
    start=time()
    # get basis names
    basisname = [name.replace('F', 'a') for name in A.columns]
        
    # perform SVD
    U, s, V = np.linalg.svd(A, full_matrices=False)  # SVD of design matrix
    W = np.diag(1/s)
    
    # get covariance matrix
    cm = np.dot(np.dot(herm(V), W**2), V)
    covmat = pd.DataFrame(cm, index=basisname, columns=basisname)
        
    # [N M] = np.shape(A);   % size of design matrix (row X column)
    
    # factors that will be used to comput regression coefficients
    Vh_dot_W = np.dot(herm(V),W)
    Uh = herm(U)
    
    try:
        instnames = b.columns
        a = {}
        for i in instnames:
            a[i] = np.dot(Vh_dot_W, np.dot(Uh, b[i]))
            
            # redchisqr = sum((A*a[i]-y[i]).^2)/(N-M);  % compute reduced chi squared
            # sa = sqrt(redchisqr*diag(Covmat));   % uncertainties in coefficients

    except AttributeError:
        # if there's only one instantiation
        a = np.dot(Vh_dot_W, np.dot(Uh, b))
        
    print time()-start
    return a, covmat
    

def results(a, cov, inv=False):
    '''
    Returns h_rec and significance for given regression parameters a and covariance
    matrix cov. If inv=True, assumes second parameter is inv(cov) instead of cov
    (for use with QR decomposition).
    '''
    h=pd.Series(index=a.keys()); s=pd.Series(index=a.keys())
    if inv==True:
        for i in a.keys():
            h.append( (abs(a[i][0]) + abs(a[i][1])) / 2 )
            s.append(np.dot(herm(a[i]), (np.dot(cov, a[i]))).real)
    else:
        for i in a.keys():
            h[i] = (abs(a[i][0]) + abs(a[i][1])) / 2
            s[i] =abs(np.dot(herm(a[i]), (np.linalg.solve(cov, a[i]))))
   
    return h, s


components = {
            'GR' : ['pl', 'cr'],
            'G4v': ['xz', 'yz'],
            'AP' : ['pl', 'cr', 'xz', 'yz', 'br'],
            'Sid': ['cos1', 'cos2', 'sin1', 'sin2', 'cnst']
            }


class Search(object):
    '''
    Searches for signals using model dependent and independent methods. Can take several
    instantiations at once (columns of pd.DF), but assumes one single PSR.
    '''
    
    def __init__(self, syst, data=False, methods=['GR', 'G4v', 'AP', 'Sid']):
        self.syst = syst
        self.det = syst.det              # LHO, LLO, etc.
        self.psr = syst.psr              # PSR name
        
        if data:
            self.initialize_data(data)
            
        self.methods = methods        

        # injection
        self.hasinjections = False
        self.haspatterns = False

        
    def initialize_data(self, data):
        self.data = data                 # pd.DF w/ index: t, columns: freqs
        self.t = np.array(data.index)
        self.instNames = data.columns
    
        # take first inst to get sigma
        if 'sg' not in dir(self):
            self.inst0 = data[self.instNames[0]]
            self.sg = 2 * getsigma(self.inst0, self.psr)
        
        # get antenna patterns only first time
        if not self.haspatterns:
            self.syst.interact(self.t, self.methods)

        # form data DF
        self.y = self.data.div(self.sg, axis=0)
        self.hinj = np.array([0]*len(self.instNames))
        
                 
    def chi(self, syst, kind):
        '''
        Chisqr regression using SVD decomposition.
        Optional: feed your own antenna patterns.
        '''
        print 'Regression to %(kind)s template...' % locals(),
        
        # in case not all
        try:
            if all([k in syst.response.kinds for k in templates.bases_aps[kind]]):
                haspatterns = True
            else:
                haspatterns = False
        except KeyError:
            haspatterns = False
            
        signal = templates.Signal(syst, kind, self.t, isloaded=haspatterns)
        signal.design_matrix()
        desMat = signal.dm
        
        A = desMat.div(self.sg, axis=0)

        a, covmat = svd_analysis(A, self.y)
        
        adf = pd.DataFrame(a, index=components[kind])
        
        print 'regression done.'

        return adf, covmat
        
        
    def qr(self, kind, F=[]):
        '''
        Chisqr regression using QR decomposition.
        Optional: feed your own antenna patterns as F=.
        '''
        print 'QR regression to %(kind)s template...' % locals(),
        
        template = simulate.Signal(self.det, 1, self.t, self.psr, F)
            
        desMat = template.design_matrix()
        A = desMat.div(self.sg, axis=0)
        
        # get covariance matrix
        invcov = A.T.dot(A)
        
        # qr decomposition of A
        Q, R = np.linalg.qr(A)
        Qt = Q.T

        # Loop over instantiations
        a={}
        try:
            for i in self.instNames:
                Qb = np.dot(Qt, self.y[i]) # computing Q^T*b (project b onto the range of A)
                a[i] = scipy.linalg.solve_triangular(R, Qb, check_finite=False) # solving R*x = Q^T*b
        except AttributeError:
            Qb = np.dot(Qt, self.y) # computing Q^T*b (project b onto the range of A)
            a = scipy.linalg.solve_triangular(R, Qb, check_finite=False) # solving R*x = Q^T*b
        
        adf = pd.DataFrame(a, index=components[kind])
        
        print 'regression done.'

        return adf, invcov

        
    def generalChi(self, data=False):
        '''
        Performs indicated searches (default: all) and returns sig and h_rec Series for
        each of them. Takes option methods = [...] indicating which regressions to do.
        
        '''
        
        if data:
            self.initialize_data(data)
            
        h = pd.DataFrame(index=self.hinj, columns=self.methods)
        s = pd.DataFrame(index=self.hinj, columns=self.methods)
        
        print 'Search methods:',
        print self.methods
        for basis in self.methods:
            a, cov = self.chi(self.syst, basis)
            h[basis], s[basis] = results(a, cov)
        
        return h, s
        
        
    def inject(self, syst, injkind, hinj, pdif):
        '''
        Injects signals in data. WARNING: original data is modified, not copied.
        Returns injection location vector.
        '''
        # get simulated signal
        try:
            if all([k in syst.response.kinds for k in templates.bases_aps[injkind]]):
                haspatterns = True
            else:
                haspatterns = False
        except KeyError:
            haspatterns = False
            
        signal = templates.Signal(syst, injkind, self.t, isloaded=haspatterns)

        signal.phases['cr'] = pdif
        signal.phases['yz'] = pdif

            
        signal.simulate(hinj)
        sim = signal.signal
        
        nfreq = len(self.data.columns)
        ninj = len(sim.columns)
        
        # build injection vector
        injLocations = [int(x) for x in np.linspace(0, nfreq, ninj, endpoint=False)]
        inj = np.zeros(nfreq)
        inj[injLocations] = hinj
        print injLocations
        # inject
        for i in range(0, ninj):
            loc  = injLocations[i]
            self.data[self.instNames[loc]] += sim[i]

        # update y vector for regressions  
        self.y = self.data.div(self.sg, axis=0)
        self.hasinjections = True
        self.hinj = inj
        
        print '%(injkind)s signal injected.' % locals()
        
