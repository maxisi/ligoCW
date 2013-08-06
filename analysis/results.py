import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import paths
from analysis import psrplot

from sys import exit
import os

reload(psrplot)

class Open(object):
    
    def __init__(self, detector, injkind, pinj, location=paths.results):
        self.detector = detector
        self.injkind = injkind
        self.pinj = pinj
        self.location = location + self.detector + '/' + self.detector + '_' + self.injkind + '_' + paths.pname(self.pinj) + '_'
        
        # load data
        self.hinj = np.load(self.location + 'hinj.npy')

        try:
            hrec_file = pd.HDFStore(self.location + 'h', 'r')
            self.hrec = {psr.strip('/'): hrec_file[psr.strip('/')] for psr in hrec_file.keys()}
        except IOError:
            print 'Cannot find hrec file in %(self.location)s' % locals()
        finally:
            hrec_file.close()
            
        try:
            s_file = pd.HDFStore(self.location + 's', 'r')
            self.s = {psr.strip('/'): s_file[psr] for psr in s_file.keys()}
        except IOError:
            print 'Cannot find s file in %(self.location)s' % locals()
        finally:
            s_file.close()
            
        self.psrlist = self.hrec.keys()
        
        
        
    def plots(self, psr, kind, method=[], extra_name=''):
    
        if method==[]:
            method=list(self.hrec[psr].columns)
            print method
        else:
            pass
        
        getattr(psrplot, kind)(hinj=self.hinj, hrec=self.hrec[psr], s=self.s[psr], methods=method)
        
        header = self.injkind + paths.pname(self.pinj) + ' injection on ' + self.detector + ' data for ' + psr + ' ' + extra_name
        plt.title(header)
        
        pltdir = 'plots/' + self.detector + '/' + self.injkind + '/' + kind + '/'
        pltname = self.detector + '_' + self.injkind + '_' + paths.pname(self.pinj) + '_' + kind + extra_name
        save_to = pltdir + pltname
        
        try:
            os.makedirs(pltdir)
        except OSError:
            pass
            
        plt.savefig(save_to, bbox_inches='tight')
        plt.close()
        
        print 'Plot saved to:\n %(save_to)s' % locals()