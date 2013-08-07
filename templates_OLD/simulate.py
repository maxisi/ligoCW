import numpy as np
import pandas as pd
from sys import exit
from templates import antennapatterns as ap
from templates import source
from templates import sidereal as sd

from time import time

components = {
            'GR' : ['F+', 'Fx'],
            'G4v': ['Fxz', 'Fyz'],
            'AP' : ['pl', 'cr', 'xz', 'yz', 'br'],
            'Sid': ['cos1', 'cos2', 'sin1', 'sin2', 'cnst']
            }

class Signal(object):
    '''
    Can generate a signal. Possible returns: signal (time series), design matrix.
    If you wish, provide your own antenna patterns setting F = {dict from tm.getAP}
    '''

    def __init__(self, det, h0, t, psr, F=[]):
        self.h0 = h0
        self.det = det
        self.psr = psr
        self.t = t
        
        # print 'Getting source information.'
        src = source.Source()
        self.inc = src.param['INC'] # get source inclinations
        self.psrlist = [psr]
        
        # print 'Getting antenna patterns information:',
        
        self.F = F
        
        self.ps = {
                    'GR' : ['pl', 'cr'],
                    'G4v': ['xz', 'yz'],
                    'AP' : ['pl', 'cr', 'xz', 'yz', 'br']
                }
        
        
    def getcomponents(self, kind):
        
        if kind not in ['GR', 'G4v', 'AP']:
            print 'ERROR: can only get componentes for GR, G4v or AP. You ask for: %r' % kind
            exit()
        else:
            pass
        
        # Relative strength catalogue
        hs = {
            'GR' : {
                    'pl': (1 + np.cos(self.inc)**2)/2.0,
                    'cr': np.cos(self.inc)
                    },
                    
            'G4v': {
                    'xz': np.sin(self.inc),
                    'yz': np.sin(self.inc)*np.cos(self.inc)
                    },
            'AP' : []
            }
            
        # Put components together
        h = pd.DataFrame(hs[kind])
        try:
        
            if all([p in self.F.keys() for p in self.ps[kind]]):
                # print 'all polarizations present'
                F = self.F
                
            else:
                # print 'some polarization was missing'
                F = ap.getAP(self.det, self.t, self.ps[kind])
                
        except AttributeError:
            # print "self.F didn't contain any polarizations"
            F = ap.getAP(self.det, self.t, self.ps[kind])
            
        return h, F
     

    def design_matrix(self, kind):
        # Returns design necessary for Chisqr regressions. Uses getcomponents to retrieve aps.

        if kind in ['GR', 'G4v']:
            h, F = self.getcomponents(kind) # gets basis and weighs (as DF) for regression
        
            i0 = self.ps[kind][0] # gets corresponding component names
            i1 = self.ps[kind][1]
            
            dm1 = F[i0].T[self.psr] * h[i0][self.psr] /2.0
            dm2 = F[i1].T[self.psr] * h[i1][self.psr] /2.0
            dm = pd.concat([dm1, dm2], axis=1, keys=components[kind]) # DF cols: comp. names, index: t.
                        
        elif kind == 'AP':
            _, F = self.getcomponents(kind)
            
            dm = pd.concat([F[p].ix[self.psr] for p in self.ps[kind]], axis=1, keys=components[kind])
            
        else: # Sid
            th = pd.Series([sd.w * int(gpst) for gpst in self.t], index=self.t)
            basis = [np.cos(th), np.cos(2*th), np.sin(th), np.sin(2*th)]
            dm = pd.concat(basis, axis=1, keys=components[kind][0:4])
            dm[components[kind][4]]=1
            
        return dm

            
    def sim(self, kind, phase=0, phasedif=np.pi/2.):
    
        i0 = self.ps[kind][0] # gets corresponding component names
        i1 = self.ps[kind][1]
        
        h, F = self.getcomponents(kind)
        dm1 = F[i0].T[self.psr] * h[i0][self.psr]
        dm2 = F[i1].T[self.psr] * h[i1][self.psr]
        
        hinj = [h_i/2.0 for h_i in self.h0]
        
        s = np.multiply.outer((dm1*np.exp(1j*phase) + dm2*np.exp(1j*(phase+phasedif))), hinj)
        
        signal = pd.DataFrame(s, index=self.t)

        return signal