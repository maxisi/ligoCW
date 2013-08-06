import numpy as np
import pandas as pd
import sidereal as sd
from time import time
from sys import exit

import paths
reload(paths)
reload(sd)

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

            
class Det(object):
    
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