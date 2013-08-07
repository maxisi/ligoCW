from time import time
from sys import exit
import matplotlib.pyplot as plt

import math
import numpy as np
import pandas as pd
import scipy.linalg

import paths
from templates_OLD import sidereal as sd
from templates_OLD import simulate
from templates import getpsrlist

from time import time


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


ps = {
        'GR' : ['pl', 'cr'],
        'G4v': ['xz', 'yz'],
        'AP' : ['pl', 'cr', 'xz', 'yz', 'br'],
        'Sid': []
    }  

components = {
            'GR' : ['GR_pl', 'GGR_cr'],
            'G4v': ['G4v_xz', 'G4v_yz'],
            'AP' : ['pl', 'cr', 'xz', 'yz', 'br'],
            'Sid': ['cos1', 'cos2', 'sin1', 'sin2', 'cnst']
            }


class Search(object):
    '''
    Searches for signals using model dependent and independent methods. Can take several
    instantiations at once (columns of pd.DF), but assumes one single PSR.
    '''
    
    def __init__(self, det, psr, data):
        self.det = det              # LHO, LLO, etc.
        self.psr = psr              # PSR name
        self.data = data            # pd.DF w/ index: t, columns: freqs
        
        self.t = np.array(data.index)
        self.instNames = data.columns
        
        # take first inst to get sigma
        self.inst0 = data[self.instNames[0]]
        self.sg = 2 * getsigma(self.inst0, self.psr)

        # form data DF
        self.y = self.data.div(self.sg, axis=0)
        
        # injection
        self.hasinjections = False
        self.hinj = np.array([0]*len(self.instNames))
        
        
    def chi(self, kind, F=[]):
        '''
        Chisqr regression using SVD decomposition.
        Optional: feed your own antenna patterns.
        '''
        print 'Regression to %(kind)s template...' % locals(),
        template = simulate.Signal(self.det, 1, self.t, self.psr, F)

        desMat = template.design_matrix(kind)

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
            
        desMat = template.design_matrix(kind)
        A = desMat.div(self.sg, axis=0)
        start = time()
        
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
        
        print time()-start
        print 'regression done.'

        return adf, invcov

        
    def generalChi(self, methods = ['Sid', 'AP', 'GR', 'G4v']):
        '''
        Performs indicated searches (default: all) and returns sig and h_rec Series for
        each of them. Takes option methods = [...] indicating which regressions to do.
        
        '''
        pol = sum([ps[basis] for basis in methods], [])
        
        if pol!=[]:
            aPs = ap.getAP(self.det, self.t, pol)
        else:
            aPs = []
            
        h = pd.DataFrame(index=self.hinj, columns=methods)
        s = pd.DataFrame(index=self.hinj, columns=methods)
             
        for basis in methods:
            a, cov = self.chi(basis, F=aPs)
            h[basis], s[basis] = results(a, cov)
        
        return h, s
        
        
    def inject(self, injkind, hinj, pdif):
        '''
        Injects signals in data. WARNING: original data is modified, not copied.
        Returns injection location vector.
        '''
        # get simulated signal
        s = simulate.Signal(self.det, hinj, self.t, self.psr)
        sim = s.sim(injkind, phasedif=pdif)
        
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
        
