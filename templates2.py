import numpy as np
import pandas as pd
from bs4 import BeautifulSoup
from urllib2 import urlopen
from collections import namedtuple

import sys, os
from sys import argv, exit

from templates import sidereal as sd
import paths


## SOURCE

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
            
            

## DETECTOR

def detnames(d):
    if d in ['H1', 'H2']:
        det = 'LHO'
    elif d == 'L1':
        det = 'LHO'
    elif d == 'V1':
        det = 'VIR'
    elif d in['LHO', 'LLO', 'VIR']:
        det = d
    else:
        print 'ERROR: %r is an unknown detector name.' % d
        exit()
    return det


# Detector parameters, all angles in radians. (Source: PRD 58, 063001 p3.)
detectors = pd.DataFrame({
        'LHO': {
                'lat': 0.8107054375513661,
                'lon': -2.084097659806429,
                'x_east': 2.199114857512855,
                'arm_ang': np.pi/2
                },
    
        'LLO': {
                'lat': 0.5333726194094671, 
                'lon': -1.584235362035253,
                'x_east': 3.4505159311927893,
                'arm_ang': np.pi/2
                },
    
        'GEO': {
                'lat': 0.9119345341670372,
                'lon': 0.17121679962064373,
                'x_east': 0.37716565135597474,
                'arm_ang': 1.646369083406251
                },
    
        'VIR': {
                'lat': 0.761487152645126,
                'lon': 0.1832595714594046,
                'x_east': 1.2479104151759457,
                'arm_ang': 1.5707963267948966
                },
    
        'TAM': {
                'lat': 0.6227334771115768,
                'lon': 2.4354324382328874,
                'x_east': 3.141592653589793,
                'arm_ang': 1.5707963267948966
                }
        })

            
class Detector(object):
    
    def __init__(self, d, t=[]):
        self.id = d  
        self.name = detnames(d)
        self.param = detectors[self.name]
        self.t = np.array([int(x) for x in t])
        self.nentries = len(t)
        self.path = paths.vectors + 'detVec_' + self.name
        
        if len(self.t)!=0:
            self.loadVectors()
        

    def fileload(self):
        '''
        Loads detector arm vectors from file.
        '''
        try:
            file = pd.HDFStore(self.path, 'r')
            self.dx = file['dx']
            self.dy = file['dy']
            self.dz = file['dz']
            file.close() 
        except IOError:
            self.createVectors()


    def loadVectors(self):
        '''
        Loads detector arm vectors if necessary.
        '''
        # Check if vectors are stored in file'
        self.fileload()
            
        # Check data type'
        try:
            self.dx.columns
        except AttributeError:
            self.fileload()
            
        if all(self.dx.index==self.t):
            # All times present.
            pass
        else:
            self.createVectors()
        

    def createVectors(self):
        '''
        Returns arm vectors in Cartesian sidereal coordinates.
        '''
        northPole = np.array([0, 0, 1])     # Earth center to North pole
        lat = self.param.ix['lat']
        lon = self.param.ix['lon']
        x_east = self.param.ix['x_east']
        arm_ang = self.param.ix['arm_ang']
        t = self.t
        length = self.nentries

        # Angle between detector and Aries (vernal equinox) at time t
        # fiducial GPS time t0=630763213 (12hUT1 1/1/2000, JD245154).
        # See http://aa.usno.navy.mil/faq/docs/GAST.php
        offset = 67310.5484088*sd.w   # Aries-Greenwich angle at fiducial time (GMST)
        th = np.add.outer(offset+sd.w*(t-630763213), lon) # (LMST) rows: t, columns: det
        
        zenith = [np.cos(lat)*np.cos(th), np.cos(lat)*np.sin(th), np.tile(np.sin(lat),(length,1))]  # [[x0, ...], [y0, ...], [z0, ...]]
        zenith /= np.sqrt(np.sum(np.array(zenith) ** 2., axis=0))
        
        localEast = np.cross(northPole,zenith, axisb=0)    # [[x, y, z], ...]
#         localEast /= np.sqrt(np.sum(localEast ** 2, axis=1))[..., None]
        
        localNorth = np.cross(zenith, localEast, axisa=0)   # [[x, y, z], ...]
        
        xArm = np.cos(x_east)*localEast+np.sin(x_east)*localNorth
        xArm /= np.sqrt(np.sum(xArm ** 2., axis=1))[..., None]
        self.dx = pd.DataFrame(xArm, index=t, columns= ['x', 'y', 'z'])
        
        perp_xz = np.cross(zenith, self.dx, axisa=0)
        yArm = xArm*np.cos(arm_ang) + perp_xz*np.sin(arm_ang) # equals perp_xz when angle between arms is 90deg
        yArm /= np.sqrt(np.sum(yArm ** 2, axis=1))[..., None]
        self.dy = pd.DataFrame(yArm, index=t, columns= ['x', 'y', 'z'])
        
        dzT = pd.DataFrame(zenith, columns=t, index= ['x', 'y', 'z'])
        self.dz = sd.rE * dzT.T
        
        try:
            f = pd.HDFStore(self.path, 'w')
            f['dx'] = self.dx
            f['dy'] = self.dy
            f['dz'] = self.dz
        finally:
            f.close()
            
            

## ANTENNA PATTERNS
class System(object):
    def __init__(self, detector, psr):
        self.detector = detector
        self.det = detnames(detector)
        self.psr = psr


bases_names = ['pl', 'cr', 'xz', 'yz', 'br']

bases_tempNames = ['GR', 'G4v', 'AP', 'Sid']

bases_aps = {
        'GR' : ['pl', 'cr'],
        'G4v': ['xz', 'yz'],
        'AP' : ['pl', 'cr', 'xz', 'yz', 'br'],
        'Sid': []
    } 
    
bases_polNames = {
            'pl': 'plus',
            'cr': 'cross',
            'xz': 'vector_x',
            'yz': 'vector_y',
            'br': 'breathing',
            'lo': 'longitudinal'
            }

# create vector container to make it easy to move vectors around
Vectors = namedtuple('Vectors', ['dx', 'dy', 'wx', 'wy', 'wz'])

# tuples indicating which vectors need to be multiplied together
# note that detector vectors should be listed second for broadcasting reasons
polComponents = {
                'pl' : [('wx','dx'), ('wx','dy'), ('wy','dx'), ('wy','dy')],
                'cr' : [('wx','dx'), ('wx','dy'), ('wy','dx'), ('wy','dy')],
                'xz' : [('wx','dx'), ('wz','dx'), ('wx','dy'), ('wz','dy')],
                'yz' : [('wy','dx'), ('wz','dx'), ('wy','dy'), ('wz','dy')],
                'br' : [('wx','dx'), ('wx','dy'), ('wy','dx'), ('wy','dy')],
                'lo' : [('wz','dx'), ('wz','dy')]
                }


class Polarizations(object):
    '''
    Produces detector response for different polarization given input detector and
    source vectors. Vectors must be in namedtuple form as defined above.
    Effectively computes dyadic product between wave and detector tensors.
    '''
    
    def __init__(self, vectors):
        # assuming vectors is of the Vectors container kind
        self.vec = vectors
    
                    
    def product(self, polKey):
        # check all necessary vector products exist. Otherwise, create them.
        for pair in polComponents[polKey]:
        
            pairName = pair[0] + pair[1]    # 'wxdx' = ('wx','dx')[0] + ('wx','dx')[1]
            
            # check if product already defined. vars(self).keys() is list of variables
            if pairName not in vars(self).keys():
            
                # get vectors
                v0 = getattr(self.vec, pair[0])     # v0 = self.vec.wx
                v1 = getattr(self.vec, pair[1])     # v1 = self.vec.dx
                
                # dot product (mind the order! detector must always be v1 to broadcast)
                setattr(self, pairName, v1.mul(v0, axis='columns').sum(axis=1))
    
    # tensor
    def plus(self):
        self.product('pl')
        pl = (self.wxdx**2 - self.wxdy**2 - self.wydx**2 + self.wydy**2)/2.
        return pl
        
    def cross(self):
        self.product('cr')
        cr = self.wxdx*self.wydx - self.wxdy*self.wydy
        return cr
    
    # vector
    def vector_x(self):
        self.product('xz')
        xz = self.wxdx*self.wzdx - self.wxdy*self.wzdy
        return xz
    
    def vector_y(self):
        self.product('yz')
        yz = self.wydx*self.wzdx - self.wydy*self.wzdy
        return yz
    
    # scalar
    def breathing(self):
        self.product('br')
        br = np.sqrt(2)*(self.wxdx**2 - self.wxdy**2 + self.wydx**2 - wydy**2)/2.
        return br
        # Added factor of sqrt(2) to distribute power equally among polarizations.
        # Same for longitudinal.
        
    def longitudinal(self):
        self.product('lo')
        lo = (np.sqrt(2)*(self.wzdx**2 - self.wzdy**2))/2.
        return lo
        # Modified:1/2 (based on derivation of dyadic products using tensors shown in 
        # "Gravitational wave polarizations" by Bryant Garcia.
        # The factor of 2 shouldn't be there)
        
    
class Response(object):
    '''
    Contains response of 'det' to signals from source 'psr' of kind 'kinds', over 't'.
    
    For default polarization angle, use 'get' method. This will try to recover patterns
    from file and will generate them if not found. It rotates source vectors by
    polarization angles established on file. After 'get' is called, the patterns are
    stored in the class, so there is no need to call again in same session.
    
    To input polarization angle, call 'create' method directly with argument 'psi='.
    This will ignore existing APs (if any) and will fall back on detector vectors to
    compute response (after rotation by 'psi'). Call with 'savefile=False' to prevent
    method from storing results. This routine will compute dyadic products once, so it
    can be called multiple times a session.
    
    To do: Add option to select psi randomly from range?
    '''
    
    def __init__(self, psr, det, t, kinds):
        self.det = detnames(det)
        self.t = np.array(t)
        
        self.psr = psr
        self.path = paths.ap + 'ap' + psr + '_' + self.det
        
        if kinds in bases_tempNames:
            # requesting bases for preset template
            self.kinds = bases_aps[kinds]
        elif all([k in basesnames for k in kinds]):
            # requesting bases by name
            if iskind(kinds, basestring):
                self.kinds = [kinds] # in case input is a single polarization
            else:
                self.kinds = kinds
        else:
            # wrong input
            print 'ERROR: %(kinds)s is not recognized as a valid basis. ap19' % locals()
            print 'Valid bases are:'
            print '\t %(basesnames)r' % locals()
            print '\t %(bases_aps)r' % locals()
            exit()
            
        self.hasvectors = False
                    
                
    def create(self, savefile=True, psi=[]):
        # Creates antenna patterns for detector detname and pulsars in list. Saves all pols.
    
        # Retrieve detector vectors
        if not self.hasvectors:
            self.src = Source(self.psr)
            self.det = Detector(self.det, self.t)
            self.hasvectors = True
            
        # Rotate source vectors
        if len(psi)==0:
            self.psi = self.src.param['POL']
        else:
            self.psi = psi
            
        wxRot = -self.src.wy*np.cos(self.psi) + self.src.wx*np.sin(self.psi)
        wyRot = self.src.wx*np.cos(self.psi) + self.src.wy*np.sin(self.psi) 
        
        # Package vectors
        vecs = Vectors(self.det.dx, self.det.dy, self.src.wx, self.src.wy, self.src.wz)      
           
        # Get polarizations
        pols = Polarizations(vecs)
        [setattr(self, k, getattr(pols, bases_polNames[k])() ) for k in self.kinds]
            
        # Save if requested
        if savefile:  
            try:
                apF = pd.HDFStore(self.path, 'w')
                for k in self.kinds:
                    apF[k] = getattr(self, k)
            finally:
                apF.close()
    

    def get(self):
        # Assumes no APs loaded. Otherwise, will re-write.
        
        # check disk
        try:
            apFile = pd.HDFStore(self.path, 'r')
            
        except IOError:
            # no patterns on file
            self.create()
            
        else:
            # file found
            try:
                filePols = [s.strip('/') for s in apFile.keys()]
            
                # Check file is alright:
                allPolsPresent = all([k in filePols for k in self.kinds])
            
                timeCoincides = set(self.t.astype(int)) == set(apFile[filePols[0]].index)
                # if file is empty, raises IndexError
            
                if allPolsPresent and timeCoincides:
                    [setattr(self, p, apFile[p]) for p in self.kinds]
                else:
                    apFile.close()
                    self.create()
                    
            except IndexError:
                apFile.close()
                self.create()
            
            finally:
                apFile.close()

    
#     def exportAPmatlab(psr, detname, t):
#         p = ap.getAP('LHO', t)
#     
#     #     psrdict = {
#     #                 'pl':pl.T[psr].tolist(),
#     #                 'cr':cr.T[psr].tolist(),
#     #                 'br':br.T[psr].tolist(),
#     #                 'lo':lo.T[psr].tolist(),
#     #                 'xz':xz.T[psr].tolist(),
#     #                 'yz':yz.T[psr].tolist()
#     #                 }
#     
#         sio.savemat('%(paths.ap)scrab_py' % locals(), p)
