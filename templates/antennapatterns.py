from templates import source, detector
import pandas as pd
import numpy as np

from sys import exit
from collections import namedtuple

import paths
import bases


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
        self.det = detector.detnames(det)
        self.t = np.array(t)
        
        self.psr = psr
        self.path = paths.ap + 'ap' + psr + '_' + self.det
        
        if kinds in bases.tempNames:
            # requesting bases for preset template
            self.kinds = bases.aps[kinds]
        elif all([k in bases.names for k in kinds]):
            # requesting bases by name
            if iskind(kinds, basestring):
                self.kinds = [kinds] # in case input is a single polarization
            else:
                self.kinds = kinds
        else:
            # wrong input
            print 'ERROR: %(kinds)s is not recognized as a valid basis. ap19' % locals()
            print 'Valid bases are:'
            print '\t %(bases.names)r' % locals()
            print '\t %(bases.aps)r' % locals()
            exit()
            
        self.hasvectors = False
                    
                
    def create(self, savefile=True, psi=[]):
        # Creates antenna patterns for detector detname and pulsars in list. Saves all pols.
    
        # Retrieve detector vectors
        if not self.hasvectors:
            self.src = source.Source(self.psr)
            self.det = detector.Det(self.det, self.t)
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
        [setattr(self, k, getattr(pols, bases.polNames[k])() ) for k in self.kinds]
            
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
