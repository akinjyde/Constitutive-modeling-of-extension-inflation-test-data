# -*- coding: utf-8 -*-
''' 
Interpolates diameter values at given pressures
Computes structural stiffness i.e. change in pressure/change in diameter
(at fixed vessel length) or its inverse i.e. compliance.

N.B. that this qty is linear and has no theoretical basis, as can be directly
inferred from exp. data
Ref (compliance): https://link.springer.com/article/10.1007%2Fs10237-014-0556-x 
'''

# Import needed standard python modules
import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.ticker as tk

def main(P_iv, OD_iv):
    
    Pa2Hg = 0.00750062
    
    # Interpolate the in vivo outer diameter at 0, 5, 10, 15, 20, and 25 mmHg
#    OD_0 = np.interp(0.0, np.sort(np.ravel(P_iv*Pa2Hg), axis=None),
#                     np.sort(np.ravel(OD_iv), axis=None))
    
    OD_2_5 = np.interp(2.5, np.sort(np.ravel(P_iv*Pa2Hg), axis=None),
                     np.sort(np.ravel(OD_iv), axis=None))
    
    OD_5 = np.interp(5.0, np.sort(np.ravel(P_iv*Pa2Hg), axis=None),
                     np.sort(np.ravel(OD_iv), axis=None))
    
    OD_7_5 = np.interp(7.5, np.sort(np.ravel(P_iv*Pa2Hg), axis=None),
                     np.sort(np.ravel(OD_iv), axis=None))
    
    OD_10 = np.interp(10.0, np.sort(np.ravel(P_iv*Pa2Hg), axis=None),
                     np.sort(np.ravel(OD_iv), axis=None))
    
    OD_15 = np.interp(15.0, np.sort(np.ravel(P_iv*Pa2Hg), axis=None),
                     np.sort(np.ravel(OD_iv), axis=None))
    
    OD_20 = np.interp(20.0, np.sort(np.ravel(P_iv*Pa2Hg), axis=None),
                     np.sort(np.ravel(OD_iv), axis=None))
    
#    OD_25 = np.interp(25.0, np.sort(np.ravel(P_iv*Pa2Hg), axis=None),
#                     np.sort(np.ravel(OD_iv), axis=None))
    
    #   Compute compliance based on Wagenseil approach i.e. secant line i.e.
    #   an average rate of change between 2 point on a curve
    
    #   Unit of the ffg is um/mmHg
    
#    C1 = (OD_2_5 - OD_0)/(2.5 - 0.0)
    C1 = (OD_5 - OD_2_5)/(5.0 - 2.5)
    C2 = (OD_7_5 - OD_5)/(7.5 - 5.0)
    C3 = (OD_10 - OD_7_5)/(10 - 7.5)
    C4 = (OD_15 - OD_10)/(15.0 - 10.0)
    C5 = (OD_20 - OD_15)/(20.0 - 15.0)
#    C7 = (OD_25 - OD_20)/(25.0 - 20.0)
    
    return (OD_2_5, OD_5, OD_7_5, OD_10, OD_15, OD_20,
            C1, C2, C3, C4, C5)
    
    
    
    
    
    
    