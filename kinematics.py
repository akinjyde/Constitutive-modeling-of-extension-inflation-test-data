# -*- coding: utf-8 -*-
''' This module computes kinematic quantities '''

# Import needed modules
from spyder_kernels.utils import iofuncs
#from spyder.utils.iofuncs import get_matlab_value
import scipy.io as sio
import numpy as np


def main(file):
    #   Read in and format provided Matlab structure array
#    data = sio.loadmat('{}.mat'.format(filename))
    data = sio.loadmat(file, appendmat=True)
    for key, value in list(data.items()):
        data[key] = iofuncs.get_matlab_value(value)
    passive = data['passive']
#    biaxial_Fl = passive['biaxial_Fl']
    biaxial_Pd = passive['biaxial_Pd']
    unloaded = passive['unloaded']
    # =========================================================================
    #   Read in unloaded geometry
    OR = 1.0e-3*unloaded['OR_mm']
    IR = 1.0e-3*unloaded['IR_mm']
    H = OR - IR
    # =========================================================================
    #   Read in and format data for the in vivo test
    test2 = biaxial_Pd['test2']
    test2data = biaxial_Pd['test2data']
    OD_iv = test2[:, 1]                      # Outer diamter [um]
    b = test2[:, 2].size
    ftiv = 1.0e-3*test2[:, 2].reshape(b, 1)  # Trans meas force [mN] to [N]
    Piv = test2data['pressurePa']            # Intraluminal pressure [Pa]
    Lzziv = test2data['lambda_z']            # Axial stretch
    iriv = 1.0e-3*test2data['ir_mm']         # Inner radius [mm] to [m]
    oriv = 1.0e-3*test2data['or_mm']         # Outer radius [mm] to [m]
    
    hiv = oriv - iriv                        # Thickness [m]
    Lqqiv = (iriv + hiv/2.0)/(IR + H/2.0)    # Circumferential stretch
    Lrriv = 1.0/(Lzziv*Lqqiv)                # Radial stretch
    Tqqiv = (Piv*iriv)/hiv                   # Circ Cauchy stress [Pa]
    Trriv = (-Piv*iriv)/(oriv + iriv)        # Radial Cauchy stress [Pa]
    fpiv = np.pi*(iriv**2)*Piv               # Force due to pressure [N]
    fiv = fpiv + ftiv                        # Total axial force [N]
    Tzziv = fiv/(np.pi*(oriv**2 - iriv**2))  # Axial Cauchy stress [Pa]
    
    # =========================================================================
    #   Read in and format data for the 4% above test
    test3 = biaxial_Pd['test3']
    test3data = biaxial_Pd['test3data']
    P_4a = test3[:, 0]                       # Intraluminal pressure [mmHg]
    OD_4a = test3[:, 1]                      # Outer diamter [um]
    c = test3[:, 2].size
    ft4a = 1.0e-3*test3[:, 2].reshape(c, 1)  # Trans meas force[mN] to [N]
    P4a = test3data['pressurePa']            # Intraluminal pressure [Pa]
    Lzz4a = test3data['lambda_z']            # Axial stretch
    ir4a = 1.0e-3*test3data['ir_mm']         # Inner radius [mm] to [m]
    or4a = 1.0e-3*test3data['or_mm']         # Outer radius [mm] to [m]
    
    h4a = or4a - ir4a                        # Thickness [m]
    Lqq4a = (ir4a + h4a/2)/(IR + H/2)        # Circumferential stretch
    Lrr4a = 1/(Lzz4a*Lqq4a)                  # Radial stretch
    Tqq4a = (P4a*ir4a)/h4a                   # Circ Cauchy stress [Pa]
    Trr4a = (-P4a*ir4a)/(or4a + ir4a)        # Radial Cauchy stress [Pa]
    fp4a = np.pi*ir4a**2*P4a                 # Force due to pressure [N]
    f4a = fp4a + ft4a                        # Total axial force [N]
    Tzz4a = f4a/(np.pi*(or4a**2 - ir4a**2))  # Axial Cauchy stress [Pa]
    
    # =========================================================================
    #   Concatenation of in vivo with 4% above    
#    Or = np.concatenate((oriv, or4a), axis=0)
#    Or.sort(axis=0)
#    Ir = np.concatenate((iriv, ir4a), axis=0)
#    Ir.sort(axis=0)
#    h = np.concatenate((hiv, h4a), axis=0)
#    h.sort(axis=0)
#    Lzz = np.concatenate((Lzziv, Lzz4a), axis=0)
#    Lzz.sort(axis=0)
#    Lqq = np.concatenate((Lqqiv, Lqq4a), axis=0)
#    Lqq.sort(axis=0)
#    Lrr = np.concatenate((Lrriv, Lrr4a), axis=0)
#    Lrr.sort(axis=0)
#    P_exp = np.concatenate((Piv, P4a), axis=0)/1e3 # converted to kPa to avoid overflow in np.exp 01/30/18
#    P_exp.sort(axis=0)
#    ft_exp = np.concatenate((ftiv, ft4a), axis=0)/1e3 # Converted to kN to avoid overflow in np.exp 01/30/18
#    ft_exp.sort(axis=0)
#    fp_exp = np.concatenate((fpiv, fp4a), axis=0)
#    fp_exp.sort(axis=0)
#    f_exp = np.concatenate((fiv, f4a), axis=0)
#    f_exp.sort(axis=0)
#    Tqq_exp = np.concatenate((Tqqiv, Tqq4a), axis=0)
#    Tqq_exp.sort(axis=0)
#    Tzz_exp = np.concatenate((Tzziv, Tzz4a), axis=0)
#    Tzz_exp.sort(axis=0)
    # =========================================================================
    #   Choose in vivo data only for fit
    Or = oriv
    Ir = iriv
    h = hiv
    Lzz = Lzziv
    Lqq = Lqqiv
    Lrr = Lrriv
    P_exp = Piv/1e3 # converted to kPa to avoid overflow in np.exp 01/13/18
    ft_exp = ftiv/1e3 # Converted to kN to avoid overflow in np.exp 01/13/18
    fp_exp = fpiv
    f_exp = fiv
    Tqq_exp = Tqqiv
    Tzz_exp = Tzziv
    # =========================================================================
    #   Choose 4% above data only for fit
    #Or = or4a
    #Ir = ir4a
    #h = h4a
    #Lzz = Lzz4a
    #Lqq = Lqq4a
    #Lrr = Lrr4a
    #P_exp = P4a
    #fp_exp = fp4a
    #ft_exp = ft4a
    #f_exp = f4a
    #Tqq_exp = Tqq4a
    #Tzz_exp = Tzz4a
    # =========================================================================
    #   fit to concatenated/or not data
    Crr = Lrr**2
    Cqq = Lqq**2
    Czz = Lzz**2
    P_iv = Piv
    
    return (Or, Ir, h, Crr, Cqq, Czz, P_exp, ft_exp, Lrriv, Lqqiv, Lzziv, iriv,
            oriv, hiv, ftiv, fiv, OD_iv, P_iv, Trriv, Tqqiv, Tzziv, Lrr4a,
            Lqq4a, Lzz4a, ir4a, or4a, h4a, ft4a, f4a, OD_4a, P_4a, Trr4a,
            Tqq4a, Tzz4a)

#
if __name__ == "__main__":
    main(file)
