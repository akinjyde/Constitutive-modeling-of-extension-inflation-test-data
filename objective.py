# -*- coding: utf-8 -*-
''' Computes the objective function to be minimized'''

# Import needed standard python modules
import numpy as np


def main(x, Crr, Cqq, Czz, Ir, Or, P_exp, ft_exp):
    
    # Import model-dependent module, and functions
    from theory_nh2ff import Pressure_th, TForce_th
        
    '''
    #   High loads:    Ferruzzi et al., 2013
    #   https://link.springer.com/article/10.1007%2Fs10439-013-0799-1
    '''
#    Ep = ((P_th(x, Crr, Cqq, Czz, Ir, Or) - P_exp)/P_exp)**2
#    Ef = ((ft_th(x, Crr, Cqq, Czz, Ir, Or, ft_exp) - ft_exp)/ft_exp)**2
    '''
    #   Compromise between low and high loads:   Ferruzzi et al., 2013
    #   https://link.springer.com/article/10.1007%2Fs10439-013-0799-1
    '''
    # [Pressure_th] = kPa, [TForce_th] = kN, see theory_'model'.py
    # [P_exp]= kPa, [ft_exp] = kN, see kinematics/kinematics_avg.py
    Ep = ((Pressure_th(x, Crr, Cqq, Czz, Ir, Or) - P_exp)/np.mean(P_exp))**2
    Ef = ((TForce_th(x, Crr, Cqq, Czz, Ir, Or) - ft_exp)/np.mean(ft_exp))**2
    '''
    #   Gleason et al., 2008 doi:10.1016/j.jbiomech.2008.08.012
    '''
#    Ep = np.sqrt(np.sum((P_th(x, Crr, Cqq, Czz, Ir, Or) - P_exp)**2)
#                 / np.sum(P_exp**2))
#    Ef = np.sqrt(np.sum((ft_th(x, Crr, Cqq, Czz, Ir, Or) - ft_exp)**2)
#                 / np.sum(ft_exp**2))
    obj = np.sum(Ep + Ef, axis=0)
    return obj
