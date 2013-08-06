import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as stats

from sys import exit

pltcolor = {
            'GR' : 'g',
            'G4v': 'r',
            'AP' : 'b',
            'Sid': 'c'
            }
            
default_style = '+'
            
def p(hinj=[], hrec=[], s=[], psrname='', detname='', style=default_style, methods=[]):
        
    for method in methods:
        # First Calculate the interquartile range
        #(http://comments.gmane.org/gmane.comp.python.scientific.user/19755)                                                                    
        data = np.sort(hrec)                                                                                                   
        upperQuartile = stats.scoreatpercentile(data,.75)                                                                      
        lowerQuartile = stats.scoreatpercentile(data,.25)                                                                      
        IQR = upperQuartile - lowerQuartile
    
    
        # Get ideal bin size
        #(http://en.wikipedia.org/wiki/Freedman%E2%80%93Diaconis_rule)
#         fdsize = 3.49*np.std(data)*len(data)**(-1./3.)
        fdsize = 2 * IQR * len(data)**(-1./3.)
            
        #Get number of bins
        #(http://stats.stackexchange.com/questions/798/calculating-optimal-number-of-bins-in-a-histogram-for-n-where-n-ranges-from-30)
        num_bins = int((np.amax(data) - np.amin(data))/fdsize)

        cumfreqs, lowlim, binsize, _ = stats.cumfreq(data, num_bins)
        pv = [1. - cdf/max(cumfreqs) for cdf in cumfreqs]
        bins = np.linspace(lowlim, num_bins*binsize, num_bins)

        plt.plot(bins, pv, style, color=pltcolor[method], label=method)
        
        plt.yscale('log')

    plt.title(detname + ' PSR ' + psrname)

    plt.xlabel('$h_{rec}$')
    plt.ylabel('1 - CDF (log scale)')

    plt.legend(numpoints=1)
    plt.savefig('plots/p_' + detname + '_' + psrname, bbox_inches='tight')
    
    print 'Plotted and saved in: ',
    print 'plots/p_' + detname + '_' + psrname
    plt.close()
    

# compound plots  
def hinjs(hinj=[], hrec=[], s=[], style=default_style, methods=[]):

    method_plot = plt.plot(hinj, s[methods], style)
        
    for i in range(0, len(method_plot)):
        plt.setp(method_plot[i], color=pltcolor[methods[i]], label=methods[i])

    plt.xlabel('$h_{inj}$')
    plt.ylabel('Significance')

    plt.legend(numpoints=1)
    

def hinjlins(hinj=[], hrec=[], s=[], style=default_style, methods=[]):

    method_plot = plt.plot(hinj, np.sqrt(s[methods]), style)
        
    for i in range(0, len(method_plot)):
        plt.setp(method_plot[i], color=pltcolor[methods[i]], label=methods[i])
    
    linsmax = max([np.amax(np.sqrt(s[m])) for m in methods])

    plt.xlabel('$h_{inj}$')
    plt.ylabel('Linearized significance')
    plt.ylim(ymax=linsmax)
    
    plt.legend(numpoints=1, loc=2)


def hinjrec(hinj=[], hrec=[], s=[], style=default_style, methods=[]):

    method_plot = plt.plot(hinj, hrec[methods], style)
    
    for i in range(0, len(method_plot)):
        plt.setp(method_plot[i], color=pltcolor[methods[i]], label=methods[i])

    hmax = max([np.amax(hrec[m]) for m in methods])
    
    plt.xlabel('$h_{inj}$')
    plt.ylabel('$h_{rec}$')
    plt.ylim(ymax=hmax)
    

# simple plots
def inj(hinj=[], hrec=[], s=[], style=default_style, methods=[]):
    plt.plot(hinj, style)
    plt.ylabel('$h_{inj}$')
    plt.xlabel('Instantiation')
    

def rec(hinj=[], hrec=[], s=[], style=default_style, methods=[]):

    method_plot = plt.plot(h[methods], style)
    
    for i in range(0, len(method_plot)):
        plt.setp(method_plot[i], color=pltcolor[methods[i]], label=methods[i])

    plt.ylabel('$h_{rec}$')
    plt.ylabel('Instantiations')

    plt.legend(numpoints=1)
    

def sig(hinj=[], hrec=[], s=[], style=default_style, methods=[]):
    
    method_plot = plt.plot(s[methods], style)
    
    for i in range(0, len(method_plot)):
        plt.setp(method_plot[i], color=pltcolor[methods[i]], label=methods[i])

    plt.ylabel('$Significance$')
    plt.ylabel('Instantiations')

    plt.legend(numpoints=1)
    

def original(detector, psr, location='files/data/'):
    try:
        d = pd.HDFStore(location + 'dataPitkin_' + detector + '.hdf5', 'r')
        a = d[psr].tolist()
        b = [x.real for x in a]
        plt.plot(b[:1000])
    finally:
        d.close()
   
    
def sigma(detector, psrname, location='files/analysis/psrSegmentSigma', compare=False):

    if compare:
        print 'comparing'
        sgS5 = pd.read_table('files/data/source/S5/sigmaS5_' + detector, names=[None, 'matlab'], sep='\s+', header=None, index_col=0)
        sgS5.plot(style='r+')
    else:
        pass
        
#     sg1 = pd.HDFStore(location + '_mat', 'r')
    
    try:
        sg = pd.HDFStore(location, 'r') #location, 'r')
        sg[psrname].plot(style='g+')
    
        plt.legend(numpoints=1)
        plt.title(psrname + detector + 'daily std')
        plt.xlabel('GPS time')
        plt.ylabel('Standard deviation')

        plt.savefig('plots/std_' + detector + '_' + psrname, bbox_inches='tight')
        print 'Plotted and saved in: ',
        print 'plots/' + detector + '_' + psrname
        plt.close()
    
    finally:
        sg.close()
                
        
def p_original(detector, psr, location='files/remote/source/'):
    d = pd.HDFStore(location + 'dataPitkin_' + detector + '.hdf5', 'r')
    a = d[psr].tolist()
    b = [abs(x) for x in a]

    # First Calculate the interquartile range
    #(http://comments.gmane.org/gmane.comp.python.scientific.user/19755)                                                                    
    data = np.sort(d[psr].tolist())                                                                                                   
    upperQuartile = stats.scoreatpercentile(data,.75)                                                                      
    lowerQuartile = stats.scoreatpercentile(data,.25)                                                                      
    IQR = upperQuartile - lowerQuartile
    
    
        # Get ideal bin size
        #(http://en.wikipedia.org/wiki/Freedman%E2%80%93Diaconis_rule)
#         fdsize = 3.49*np.std(data)*len(data)**(-1./3.)
    fdsize = 2 * IQR * len(data)**(-1./3.)
        
    #Get number of bins
    #(http://stats.stackexchange.com/questions/798/calculating-optimal-number-of-bins-in-a-histogram-for-n-where-n-ranges-from-30)
    num_bins = int((np.amax(data) - np.amin(data))/fdsize)

    cumfreqs, lowlim, binsize, _ = stats.cumfreq(data, num_bins)
    pv = [1. - cdf/max(cumfreqs) for cdf in cumfreqs]
    bins = np.linspace(lowlim, num_bins*binsize, num_bins)

    plt.plot(bins, pv, style, color=pltcolor[method], label=method)
    
    plt.yscale('log')

    plt.title(detname + ' PSR ' + psrname)

    plt.xlabel('$h_{rec}$')
    plt.ylabel('1 - CDF (log scale)')

    plt.legend(numpoints=1)
    plt.savefig('plots/p_' + detname + '_' + psrname, bbox_inches='tight')
    
    print 'Plotted and saved in: ',
    print 'plots/p_' + detname + '_' + psrname
    plt.close()