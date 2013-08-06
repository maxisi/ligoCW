import pandas as pd
import numpy as np
import time
import sys
import os

from sys import argv, exit
from bs4 import BeautifulSoup

from urllib2 import urlopen

import sidereal as sd
import paths

paramNames = ['#', None, 'RAS', 'RAS error', 'DEC', 'DEC error']
extraParamNames = [None, 'POL', 'POL error', 'INC', 'INC error']

def currentCatalogue(display=False):
    '''
    Returns contents of current pulsar catalogue. Mainly for maintenance reasons.
    '''
    try:
        cat = pd.load(paths.psrcat)
        if display:
            print cat
        return cat
    except:
        print 'No pulsar catalogue found.'
        exit()
        
    
def getpsrlist():
    try:
        with open(paths.psrlist, 'r') as f:
            l = f.read()
        
            if len(l) == 0:
                raise IOError
            else:
                psr_list = l.split(',')
                return psr_list
            
    except IOError:
        print 'ERROR: no pulsar list!'
        exit(0)
        

class Source(object):
    '''
    Assuming single pulsar.
    '''
    
    def __init__(self, psr):
            
        self.psr = psr
        self.npsrs = 1
        self.path = paths.vectors + 'srcVec_' + psr
                     
        # If necessary, rebuild catalogue; otherwise, just load catalogue.
        try:
            f = open(paths.textfromATNF, 'r')
            pd.load(paths.psrcat)
        except IOError:
            self.build_catalogue()
            f = open(paths.textfromATNF, 'r')
        finally:
            f_text = f.read()
            f.close()
        
        if self.psr not in f_text:
            self.build_catalogue()

        psrcat = pd.load(paths.psrcat)
        self.param = psrcat.ix[self.psr]
        
        # Load vectors if a time input was given
        self.loadVectors()
        

    def build_catalogue(self, extrapsrs=[]):
        '''
        Gets location parameters for pulsars in input list from ATNF online catalogue.
        Creates and pickles corresponding pandas DataFrame 'pulsar_parameters'.
        '''
        psrs = [self.psr] + extrapsrs
        
        def atnfurl(psr_list):
            pre = 'http://www.atnf.csiro.au/research/pulsar/psrcat/proc_form.php?version=1.47&table_top.x=-438&table_top.y=-368&JName=JName&RaJ=RaJ&DecJ=DecJ&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val=&sort_attr=jname&sort_order=asc&condition=&'
            post = '&ephemeris=selected&coords_unit=raj%2Fdecj&radius=&coords_1=&coords_2=&style=Long+with+errors&no_value=*&nohead=nohead&fsize=3&x_axis=&x_scale=linear&y_axis=&y_scale=linear&state=query'
            names= 'pulsar_names='
    
            for psr in psr_list:
                names+=psr.replace('+', '%2B')
        
                if psr != psr_list[-1]:
                    names+='%0D%0A'
            url = "%(pre)s%(names)s%(post)s" % locals()
    
            if len(url)>2000:
                print 'WARNING! URL %d characters!' % len(url)
    
            return url

        # Get data
        url = atnfurl(psrs)
        soup = BeautifulSoup(urlopen(url)) # get webpage and parse
        text = str(soup.pre.string)

        # Write txt file 
        f = open(paths.textfromATNF, 'w')
        f.write(text)
        f.close()

        # Load and create DataFrame
        psr = pd.read_table(paths.textfromATNF, sep='\s+', comment='*', names= paramNames, header=None, skiprows=1, index_col=1)
        psrcat=psr.drop("#", axis=1)
    
        # Format
        formatRAS = lambda x: sd.hms_rad(x)
        formatDEC = lambda y: np.radians(sd.dms_deg(y))
        psrcat['RAS'] = psrcat['RAS'].map(formatRAS)
        psrcat['DEC'] = psrcat['DEC'].map(formatDEC)
    
        #Check extra parameters
        extra = pd.read_table(paths.psrextra, sep=',', header=None, index_col=[0], names=extraParamNames)
        psrcatComplete = pd.merge(psrcat,extra, left_index=True, right_index=True, how='outer')
    
        psrcatComplete['POL error']=psrcatComplete['POL error'].fillna(value=np.pi/4.)
        psrcatComplete['INC error']=psrcatComplete['INC error'].fillna(value=np.pi/4.)

    
        psrcatComplete.fillna(value=0, inplace=True)
        psrcatComplete.save(paths.psrcat)
     
     
    def loadVectors(self):
        '''
        Loads detector arm vectors from file.
        '''
        try:
            file = pd.HDFStore(self.path, 'r')
            self.wx = file['wx']
            self.wy = file['wy']
            self.wz = file['wz']
            file.close()
        except IOError:
            self.createVectors()
           
            
            
    def createVectors(self):
        # Return wave vectors for all sources listed.
        
        north = np.array([0, 0, 1])
        
        # take source location vector components in celestial coordinates and invert direction multiplying by -1 to get wave vector wz
        wz = [-np.cos(self.param['DEC'])*np.cos(self.param['RAS']), -np.cos(self.param['DEC'])*np.sin(self.param['RAS']), -np.sin(self.param['DEC'])] 
        self.wz = pd.Series(wz, name=self.psr, index=['x', 'y', 'z'])
        
        
        wy = np.cross(wz, north) 
        wy /= np.sqrt(np.sum(wy ** 2))
        self.wy = pd.Series(wy, name=self.psr, index=['x','y','z'])

        wx = np.cross(wy, wz)
        wx /= np.sqrt(np.sum(wx ** 2))
        self.wx = pd.Series(wx, name=self.psr, index=['x','y','z'])
                
        try:
            f = pd.HDFStore(self.path, 'w')
            f['wx'] = self.wx
            f['wy'] = self.wy
            f['wz'] = self.wz
        finally:
            f.close()