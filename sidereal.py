import numpy as np

# Seconds in a sidereal day
ss = 86164.0905

# Sidereal angular frequency of Earth
w = 2*np.pi/ss

# Earth radius (m)
rE = 6378.137e3

# LIGO data sampling period
periodLIGO = 1/16384. # in seconds, from M. Pitkin

# Speed of light (m/s)
c = 299792458.


def hmsformat(*args):
    
    if len(args) == 1:
        # Assume hh:mm:ss format
        if type(args) != str:
            argument = args[0][0]
        else:
            argument = args
            
        hms = argument.split(':')
        h, m, s = [float(x) for x in hms]
        
    elif len(args) == 3:
        h, m, s = [x for x in args]
    else:
        print 'ERROR in hmsformat: can\'t take %d arguments' % len(args)    
    return h, m, s
    

def hms_rad(*args):
    # Converts hours, minutes, seconds to radians using the sidereal angular frequency of the Earth
    h, m, s = hmsformat(args)
    sec = s + 60*(m + 60*h)
    return sec*w
    
    
def dms_deg(*args):
    # Converts degrees, minutes, seconds to decimal degrees
    d, m, s = hmsformat(args)
    return d + m/60 + s/(60**2)