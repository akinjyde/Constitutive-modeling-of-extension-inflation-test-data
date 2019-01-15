# -*- coding: utf-8 -*-
''' This module performs optimization via the Differential Evolution algorithm '''

# Import needed standard python modules
import numpy as np
from scipy.optimize import differential_evolution
# import scipy.io as sio


# Import needed created modules
from objective import main as runObj
# Import model-dependent module, and functions
from theory_nh2ff import Pressure_th, TForce_th


def main(bounds2, Crr, Cqq, Czz, Ir, Or, P_exp, ft_exp, filename, path):
    """Perform global optimization on P-d_f-l test data via the differential
    evolution function.
    
    Parameters
    ----------
    bounds2 : sequence of lower and upper bounds
    
    Crr, Cqq, Czz : diagonal elements of the right/left Cauchy-Green Strain
    tensor
    
    Ir, Or : inner, and outer radii
    
    Returns
    -------
    resDE : optmized parameters, and other details 
    """
    print('Performing Differential Evolution optmization ...')
    resDE = differential_evolution(runObj, bounds=bounds2,
                                   args=(Crr, Cqq, Czz, Ir, Or, P_exp, ft_exp),
                                   strategy='best1bin',maxiter=1000, popsize=15,
                                   tol=0.01, mutation=(0.5,1), recombination=0.9,
                                   seed=None, callback=None, disp=False,
                                   polish=True, init='latinhypercube', atol=0)
    
    print('Optimized Par [Diff_Evol]: [    mu    c1      c2      alpha]')
    
    np.set_printoptions(precision=4)    
    print('Optimized Par [Diff_Evol]:', resDE.x)
    
    print('Results for NH_2FF')
        
    # Compute the Goodness-of-fit value
    SSEp = np.sum((P_exp - Pressure_th(resDE.x, Crr, Cqq, Czz, Ir, Or))**2)
    SSTp = np.sum((P_exp - np.mean(P_exp))**2)

    SSEf = np.sum((ft_exp - TForce_th(resDE.x, Crr, Cqq, Czz, Ir, Or))**2)
    SSTf = np.sum((ft_exp - np.mean(ft_exp))**2)

    gof = np.around(1.0 - ((SSEp + SSEf)/(SSTp + SSTf)), decimals=4)
    print('GOF [Diff_Evol]:', gof)
    
#    mdict = {'resDE': resDE, 'GOF': gof}
#    sio.savemat("%s\%s\%s\%s" % (path, filename, model, 'Opt_results'), mdict=mdict, appendmat=True)

    return resDE, gof


