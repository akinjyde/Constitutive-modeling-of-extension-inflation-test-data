# -*- coding: utf-8 -*-
''' Computes the theoretical pressure and transducer measured force '''

import numpy as np
#from kinematics import main as runKin

def Pressure_th(x, Crr, Cqq, Czz, Ir, Or):
    """Computes the theoretical intraluminal pressure.
    
    Parameters
    ----------
    x : vector of model parameters
    
    Crr, Cqq, Czz : diagonal elements of the right/left Cauchy-Green Strain tensor
    
    Ir, Or : inner, and outer radii
    
    Returns
    -------
    P_th : Intraluminal pressure 
    """
    #      x[0]=ce   x[1]=c1   x[2]=c2   x[3]=alpha
    
    I4 = Cqq*(np.sin(x[3]))**2 + Czz*(np.cos(x[3]))**2
    W1 = 0.5*x[0]
    W4 = 2.0*x[1]*(I4 - 1.0)*np.exp(x[2]*(I4 - 1.0)**2)
    P_th = 2.0*(W1*(Cqq - Crr) + W4*Cqq*(np.sin(x[3]))**2)*np.log(Or/Ir)  # kPa
    return P_th


def TForce_th(x, Crr, Cqq, Czz, Ir, Or):
    """Computes the theoretical transducer-measured force.
    
    Parameters
    ----------
    x : vector of model parameters
    
    Crr, Cqq, Czz : diagonal elements of the right/left Cauchy-Green Strain tensor
    
    Ir, Or : inner, and outer radii
    
    Returns
    -------
    ft_th : Transducer-measured force 
    """
    I4 = Cqq*(np.sin(x[3]))**2 + Czz*(np.cos(x[3]))**2
    W1 = 0.5*x[0]
    W4 = 2.0*x[1]*(I4 - 1.0)*np.exp(x[2]*(I4 - 1.0)**2)
    # int => integrand
    f_int = 2.0*(W1*(2.0*Czz - Crr - Cqq)
                 + W4*(2.0*Czz*(np.cos(x[3]))**2 - Cqq*(np.sin(x[3]))**2))
    ft_th = 0.5*np.pi*(Or**2 - Ir**2)*f_int                 # kN i.e. m^2 * kPa
    return ft_th



