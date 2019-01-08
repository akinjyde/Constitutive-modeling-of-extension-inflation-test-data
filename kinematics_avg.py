# -*- coding: utf-8 -*-
''' This module computes the average kinematic quantities '''

# Import needed modules
from spyder.utils.iofuncs import get_matlab_value
import scipy.io as sio
import numpy as np
from scipy.stats import sem


def main(file1, file2, file3, file4, file5, file6, file7, file8, path, filename, model):
    
    # Read in and format provided Matlab structure array
    data1 = sio.loadmat(file1, appendmat=True)
    for key, value in list(data1.items()):
        data1[key] = get_matlab_value(value)
    passive1 = data1['passive']
#    biaxial_Fl1 = passive1['biaxial_Fl']
    biaxial_Pd1 = passive1['biaxial_Pd']
    unloaded1 = passive1['unloaded']
    # =========================================================================
    #   Read in unloaded geometry
    OR1 = 1.0e-3*unloaded1['OR_mm']             # Outer (unloaded) radius [mm]
    IR1 = 1.0e-3*unloaded1['IR_mm']             # Inner (unloaded) radius [mm]
    H1 = OR1 - IR1                              # Thickness (unloaded) [mm]
    # =========================================================================
    #   Read in and format data for the in vivo test
    test2_1 = biaxial_Pd1['test2']
    test2data1 = biaxial_Pd1['test2data']
    OD_1 = test2_1[:, 1]                        # Outer diamter [um]
    b1 = test2_1[:, 2].size
    ft1 = 1.0e-3*test2_1[:, 2].reshape(b1, 1)   # Trans meas force[mN] to [N]
    P1 = test2data1['pressurePa']               # Intraluminal pressure [Pa]
    Lzz1 = test2data1['lambda_z']               # Axial stretch
    ir1 = 1.0e-3*test2data1['ir_mm']            # Inner radius [mm] to [m]
    or1 = 1.0e-3*test2data1['or_mm']            # Outer radius [mm] to [m]
    
    h1 = or1 - ir1                              # Thickness [m]
    Lqq1 = (ir1 + h1/2.0)/(IR1 + H1/2.0)        # Circumferential stretch
    Lrr1 = 1.0/(Lzz1*Lqq1)                      # Radial stretch
    Tqq1 = (P1*ir1)/h1                          # Circ Cauchy stress [Pa]
    Trr1 = (-P1*ir1)/(or1 + ir1)                # Radial Cauchy stress [Pa]
    fp1 = np.pi*(ir1**2)*P1                     # Force due to pressure [N]
    f1 = fp1 + ft1                              # Total axial force [N]
    Tzz1 = f1/(np.pi*(or1**2 - ir1**2))         # Axial Cauchy stress [Pa]
    # =========================================================================
    #   Read in and format data for the 4% Above test
    test2_1a = biaxial_Pd1['test3']
    test2data1a = biaxial_Pd1['test3data']
    OD_1a = test2_1a[:, 1]                      # Outer diamter [um]
    b1a = test2_1a[:, 2].size
    ft1a = 1.0e-3*test2_1a[:, 2].reshape(b1a, 1)# Trans meas force[mN] to [N]
    P1a = test2data1a['pressurePa']             # Luminal pressure [Pa]
    Lzz1a = test2data1a['lambda_z']             # Axial stretch
    ir1a = 1.0e-3*test2data1a['ir_mm']          # Inner radius [mm] to [m]
    or1a = 1.0e-3*test2data1a['or_mm']          # Outer radius [mm] to [m]
    
    h1a = or1a - ir1a                           # Thickness [m]
    Lqq1a = (ir1a + h1a/2.0)/(IR1 + H1/2.0)     # Circumferential stretch
    Lrr1a = 1.0/(Lzz1a*Lqq1a)                   # Radial stretch
    Tqq1a = (P1a*ir1a)/h1a                      # Circ Cauchy stress [Pa]
    Trr1a = (-P1a*ir1a)/(or1a + ir1a)           # Radial Cauchy stress [Pa]
    fp1a = np.pi*(ir1a**2)*P1a                  # Force due to pressure [N]
    f1a = fp1a + ft1a                           # Total axial force [N]
    Tzz1a = f1a/(np.pi*(or1a**2 - ir1a**2))     # Axial Cauchy stress [Pa]
    # =========================================================================
    #   Read in and format provided Matlab structure array
    data2 = sio.loadmat(file2, appendmat=True)
    for key, value in list(data2.items()):
        data2[key] = get_matlab_value(value)
    passive2 = data2['passive']
#    biaxial_Fl2 = passive2['biaxial_Fl']
    biaxial_Pd2 = passive2['biaxial_Pd']
    unloaded2 = passive2['unloaded']
    # =========================================================================
    #   Read in unloaded geometry
    OR2 = 1.0e-3*unloaded2['OR_mm']
    IR2 = 1.0e-3*unloaded2['IR_mm']
    H2 = OR2 - IR2
    # =========================================================================
    #   Read in and format data for the in vivo test
    test2_2 = biaxial_Pd2['test2']
    test2data2 = biaxial_Pd2['test2data']
    OD_2 = test2_2[:, 1]                        # Outer diamter [um]
    b2 = test2_2[:, 2].size
    ft2 = 1.0e-3*test2_2[:, 2].reshape(b2, 1)   # Trans meas force[mN] to [N]
    P2 = test2data2['pressurePa']               # Luminal pressure [Pa]
    Lzz2 = test2data2['lambda_z']               # Axial stretch
    ir2 = 1.0e-3*test2data2['ir_mm']            # Inner radius [mm] to [m]
    or2 = 1.0e-3*test2data2['or_mm']            # Outer radius [mm] to [m]
    
    h2 = or2 - ir2                              # Thickness [m]
    Lqq2 = (ir2 + h2/2.0)/(IR2 + H2/2.0)        # Circ stretch
    Lrr2 = 1.0/(Lzz2*Lqq2)                      # Radial stretch
    Tqq2 = (P2*ir2)/h2                          # Circ Cauchy stress [Pa]
    Trr2 = (-P2*ir2)/(or2 + ir2)                # Radial Cauchy stress [Pa]
    fp2 = np.pi*(ir2**2)*P2                     # Force due to pressure [N]
    f2 = fp2 + ft2                              # Total axial force [N]
    Tzz2 = f2/(np.pi*(or2**2 - ir2**2))         # Axial Cauchy stress [Pa]
    # =========================================================================
    #   Read in and format data for the 4% ABOVE test
    test2_2a = biaxial_Pd2['test3']
    test2data2a = biaxial_Pd2['test3data']
    OD_2a = test2_2a[:, 1]                      # Outer diamter [um]
    b2a = test2_2a[:, 2].size
    ft2a = 1.0e-3*test2_2a[:, 2].reshape(b2a, 1)# Trans meas force[mN] to [N]
    P2a = test2data2a['pressurePa']             # Luminal pressure [Pa]
    Lzz2a = test2data2a['lambda_z']             # Axial stretch
    ir2a = 1.0e-3*test2data2a['ir_mm']          # Inner radius [mm] to [m]
    or2a = 1.0e-3*test2data2a['or_mm']          # Outer radius [mm] to [m]
    
    h2a = or2a - ir2a                           # Thickness [m]
    Lqq2a = (ir2a + h2a/2.0)/(IR2 + H2/2.0)     # Circumferential stretch
    Lrr2a = 1.0/(Lzz2a*Lqq2a)                   # Radial stretch
    Tqq2a = (P2a*ir2a)/h2a                      # Circ Cauchy stress [Pa]
    Trr2a = (-P2a*ir2a)/(or2a + ir2a)           # Radial Cauchy stress [Pa]
    fp2a = np.pi*(ir2a**2)*P2a                  # Force due to pressure [N]
    f2a = fp2a + ft2a                           # Total axial force [N]
    Tzz2a = f2a/(np.pi*(or2a**2 - ir2a**2))     # Axial Cauchy stress [Pa]
    # =========================================================================
    #   Read in and format provided Matlab structure array
    data3 = sio.loadmat(file3, appendmat=True)
    for key, value in list(data3.items()):
        data3[key] = get_matlab_value(value)
    passive3= data3['passive']
#    biaxial_Fl3= passive3['biaxial_Fl']
    biaxial_Pd3= passive3['biaxial_Pd']
    unloaded3= passive3['unloaded']
    # =========================================================================
    #   Read in unloaded geometry
    OR3= 1.0e-3*unloaded3['OR_mm']
    IR3= 1.0e-3*unloaded3['IR_mm']
    H3= OR3- IR3
    # =========================================================================
    #   Read in and format data for the in vivo test
    test2_3= biaxial_Pd3['test2']
    test2data3= biaxial_Pd3['test2data']
    OD_3= test2_3[:, 1]                         # Outer diamter [um]
    b3= test2_3[:, 2].size
    ft3= 1.0e-3*test2_3[:, 2].reshape(b3, 1)    # Trans meas force[mN] to [N]
    P3= test2data3['pressurePa']                # Intraluminal pressure [Pa]
    Lzz3= test2data3['lambda_z']                # Axial stretch
    ir3= 1.0e-3*test2data3['ir_mm']             # Inner radius [mm] to [m]
    or3= 1.0e-3*test2data3['or_mm']             # Outer radius [mm] to [m]
    
    h3= or3- ir3                                # Thickness [m]
    Lqq3= (ir3+ h3/2.0)/(IR3+ H3/2.0)           # Circumferential stretch
    Lrr3= 1.0/(Lzz3*Lqq3)                       # Radial stretch
    Tqq3= (P3*ir3)/h3                           # Circ Cauchy stress [Pa]
    Trr3= (-P3*ir3)/(or3+ ir3)                  # Radial Cauchy stress [Pa]
    fp3= np.pi*(ir3**2)*P3                      # Force due to pressure [N]
    f3= fp3+ ft3                                # Total axial force [N]
    Tzz3= f3/(np.pi*(or3**2 - ir3**2))          # Axial Cauchy stress [Pa]
    # =========================================================================
    #   Read in and format data for 4% ABOVE test
    test2_3a= biaxial_Pd3['test3']
    test2data3a= biaxial_Pd3['test3data']
    OD_3a= test2_3a[:, 1]                       # Outer diamter [um]
    b3a= test2_3a[:, 2].size
    ft3a= 1.0e-3*test2_3a[:, 2].reshape(b3a, 1) # Trans meas force[mN] to [N]
    P3a= test2data3a['pressurePa']              # Intraluminal pressure [Pa]
    Lzz3a= test2data3a['lambda_z']              # Axial stretch
    ir3a= 1.0e-3*test2data3a['ir_mm']           # Inner radius [mm] to [m]
    or3a= 1.0e-3*test2data3a['or_mm']           # Outer radius [mm] to [m]
    
    h3a= or3a- ir3a                             # Thickness [m]
    Lqq3a= (ir3a+ h3a/2.0)/(IR3+ H3/2.0)        # Circumferential stretch
    Lrr3a= 1.0/(Lzz3a*Lqq3a)                    # Radial stretch
    Tqq3a= (P3a*ir3a)/h3a                       # Circ Cauchy stress [Pa]
    Trr3a= (-P3a*ir3a)/(or3a+ ir3a)             # Radial Cauchy stress [Pa]
    fp3a= np.pi*(ir3a**2)*P3a                   # Force due to pressure [N]
    f3a= fp3a + ft3a                            # Total axial force [N]
    Tzz3a= f3a/(np.pi*(or3a**2 - ir3a**2))      # Axial Cauchy stress [Pa]
    # =========================================================================
    #   Read in and format provided Matlab structure array
    data4 = sio.loadmat(file4, appendmat=True)
    for key, value in list(data4.items()):
        data4[key] = get_matlab_value(value)
    passive4 = data4['passive']
#    biaxial_Fl4 = passive4['biaxial_Fl']
    biaxial_Pd4 = passive4['biaxial_Pd']
    unloaded4 = passive4['unloaded']
    # =========================================================================
    #   Read in unloaded geometry
    OR4 = 1.0e-3*unloaded4['OR_mm']
    IR4 = 1.0e-3*unloaded4['IR_mm']
    H4 = OR4- IR4
    # =========================================================================
    #   Read in and format data for the in vivo test
    test2_4 = biaxial_Pd4['test2']
    test2data4 = biaxial_Pd4['test2data']
    OD_4 = test2_4[:, 1]                        # Outer diamter [um]
    b4 = test2_4[:, 2].size
    ft4 = 1.0e-3*test2_4[:, 2].reshape(b4, 1)   # Trans meas force[mN] to [N]
    P4 = test2data4['pressurePa']               # Intraluminal pressure [Pa]
    Lzz4 = test2data4['lambda_z']               # Axial stretch
    ir4 = 1.0e-3*test2data4['ir_mm']            # Inner radius [mm] to [m]
    or4 = 1.0e-3*test2data4['or_mm']            # Outer radius [mm] to [m]
    
    h4 = or4- ir4                               # Thickness [m]
    Lqq4 = (ir4+ h4/2.0)/(IR4+ H4/2.0)          # Circumferential stretch
    Lrr4 = 1.0/(Lzz4*Lqq4)                      # Radial stretch
    Tqq4 = (P4*ir4)/h4                          # Circ Cauchy stress [Pa]
    Trr4 = (-P4*ir4)/(or4+ ir4)                 # Radial Cauchy stress [Pa]
    fp4 = np.pi*(ir4**2)*P4                     # Force due to pressure [N]
    f4 = fp4+ ft4                               # Total axial force [N]
    Tzz4 = f4/(np.pi*(or4**2 - ir4**2))         # Axial Cauchy stress [Pa]
    # =========================================================================
    #   Read in and format data for 4% ABOVE test
    test2_4a = biaxial_Pd4['test3']
    test2data4a = biaxial_Pd4['test3data']
    OD_4a = test2_4a[:, 1]                      # Outer diamter [um]
    b4a = test2_4a[:, 2].size
    ft4a = 1.0e-3*test2_4a[:, 2].reshape(b4a, 1)# Trans meas force[mN] to [N]
    P4a = test2data4a['pressurePa']             # Intraluminal pressure [Pa]
    Lzz4a = test2data4a['lambda_z']             # Axial stretch
    ir4a = 1.0e-3*test2data4a['ir_mm']          # Inner radius [mm] to [m]
    or4a = 1.0e-3*test2data4a['or_mm']          # Outer radius [mm] to [m]
    
    h4a = or4a- ir4a                            # Thickness [m]
    Lqq4a = (ir4a + h4a/2.0)/(IR4+ H4/2.0)      # Circumferential stretch
    Lrr4a = 1.0/(Lzz4a*Lqq4a)                   # Radial stretch
    Tqq4a = (P4a*ir4a)/h4a                      # Circ Cauchy stress [Pa]
    Trr4a = (-P4a*ir4a)/(or4a + ir4a)           # Radial Cauchy stress [Pa]
    fp4a = np.pi*(ir4a**2)*P4a                  # Force due to pressure [N]
    f4a = fp4a + ft4a                           # Total axial force [N]
    Tzz4a = f4a/(np.pi*(or4a**2 - ir4a**2))     # Axial Cauchy stress [Pa]
    # =========================================================================
    #   Read in and format provided Matlab structure array
    data5 = sio.loadmat(file5, appendmat=True)
    for key, value in list(data5.items()):
        data5[key] = get_matlab_value(value)
    passive5 = data5['passive']
#    biaxial_Fl5 = passive5['biaxial_Fl']
    biaxial_Pd5 = passive5['biaxial_Pd']
    unloaded5 = passive5['unloaded']
    # =========================================================================
    #   Read in unloaded geometry
    OR5 = 1.0e-3*unloaded5['OR_mm']
    IR5 = 1.0e-3*unloaded5['IR_mm']
    H5 = OR5- IR5
    # =========================================================================
    #   Read in and format data for the in vivo test
    test2_5 = biaxial_Pd5['test2']
    test2data5 = biaxial_Pd5['test2data']
    OD_5 = test2_5[:, 1]                        # Outer diamter [um]
    b5 = test2_5[:, 2].size
    ft5 = 1.0e-3*test2_5[:, 2].reshape(b5, 1)   # Trans meas force[mN] to [N]
    P5 = test2data5['pressurePa']               # Intraluminal pressure [Pa]
    Lzz5 = test2data5['lambda_z']               # Axial stretch
    ir5 = 1.0e-3*test2data5['ir_mm']            # Inner radius [mm] to [m]
    or5 = 1.0e-3*test2data5['or_mm']            # Outer radius [mm] to [m]
    
    h5 = or5 - ir5                              # Thickness [m]
    Lqq5 = (ir5 + h5/2.0)/(IR5 + H5/2.0)        # Circumferential stretch
    Lrr5 = 1.0/(Lzz5*Lqq5)                      # Radial stretch
    Tqq5 = (P5*ir5)/h5                          # Circ Cauchy stress [Pa]
    Trr5 = (-P5*ir5)/(or5 + ir5)                # Radial Cauchy stress [Pa]
    fp5 = np.pi*(ir5**2)*P5                     # Force due to pressure [N]
    f5 = fp5 + ft5                              # Total axial force [N]
    Tzz5 = f5/(np.pi*(or5**2 - ir5**2))         # Axial Cauchy stress [Pa]
    # =========================================================================
    #   Read in and format data for 4% ABOVE test
    test2_5a = biaxial_Pd5['test3']
    test2data5a = biaxial_Pd5['test3data']
    OD_5a = test2_5a[:, 1]                      # Outer diamter [um]
    b5a = test2_5a[:, 2].size
    ft5a = 1.0e-3*test2_5a[:, 2].reshape(b5a, 1)# Trans meas force[mN] to [N]
    P5a = test2data5a['pressurePa']             # Intraluminal pressure [Pa]
    Lzz5a = test2data5a['lambda_z']             # Axial stretch
    ir5a = 1.0e-3*test2data5a['ir_mm']          # Inner radius [mm] to [m]
    or5a = 1.0e-3*test2data5a['or_mm']          # Outer radius [mm] to [m]
    
    h5a = or5a - ir5a                           # Thickness [m]
    Lqq5a = (ir5a + h5a/2.0)/(IR5 + H5/2.0)     # Circumferential stretch
    Lrr5a = 1.0/(Lzz5a*Lqq5a)                   # Radial stretch
    Tqq5a = (P5a*ir5a)/h5a                      # Circ Cauchy stress [Pa]
    Trr5a = (-P5a*ir5a)/(or5a + ir5a)           # Radial Cauchy stress [Pa]
    fp5a = np.pi*(ir5a**2)*P5a                  # Force due to pressure [N]
    f5a = fp5a + ft5a                           # Total axial force [N]
    Tzz5a = f5a/(np.pi*(or5a**2 - ir5a**2))     # Axial Cauchy stress [Pa]
    # =========================================================================
    #   Read in and format provided Matlab structure array
    data6 = sio.loadmat(file6, appendmat=True)
    for key, value in list(data6.items()):
        data6[key] = get_matlab_value(value)
    passive6 = data6['passive']
#    biaxial_Fl6 = passive6['biaxial_Fl']
    biaxial_Pd6 = passive6['biaxial_Pd']
    unloaded6 = passive6['unloaded']
    # =========================================================================
    #   Read in unloaded geometry
    OR6 = 1.0e-3*unloaded6['OR_mm']
    IR6 = 1.0e-3*unloaded6['IR_mm']
    H6 = OR6- IR6
    # =========================================================================
    #   Read in and format data for the in vivo test
    test2_6 = biaxial_Pd6['test2']
    test2data6 = biaxial_Pd6['test2data']
    OD_6 = test2_6[:, 1]                        # Outer diamter [um]
    b6 = test2_6[:, 2].size
    ft6 = 1.0e-3*test2_6[:, 2].reshape(b6, 1)   # Trans meas force[mN] to [N]
    P6 = test2data6['pressurePa']               # Intraluminal pressure [Pa]
    Lzz6 = test2data6['lambda_z']               # Axial stretch
    ir6 = 1.0e-3*test2data6['ir_mm']            # Inner radius [mm] to [m]
    or6 = 1.0e-3*test2data6['or_mm']            # Outer radius [mm] to [m]
    
    h6 = or6 - ir6                              # Thickness [m]
    Lqq6 = (ir6 + h6/2.0)/(IR6 + H6/2.0)        # Circumferential stretch
    Lrr6 = 1.0/(Lzz6*Lqq6)                      # Radial stretch
    Tqq6 = (P6*ir6)/h6                          # Circ Cauchy stress [Pa]
    Trr6 = (-P6*ir6)/(or6 + ir6)                # Radial Cauchy stress [Pa]
    fp6 = np.pi*(ir6**2)*P6                     # Force due to pressure [N]
    f6 = fp6 + ft6                              # Total axial force [N]
    Tzz6 = f6/(np.pi*(or6**2 - ir6**2))         # Axial Cauchy stress [Pa]
    # =========================================================================
    #   Read in and format data for 4% ABOVE test
    test2_6a = biaxial_Pd6['test3']
    test2data6a = biaxial_Pd6['test3data']
    OD_6a = test2_6a[:, 1]                      # Outer diamter [um]
    b6a = test2_6a[:, 2].size
    ft6a = 1.0e-3*test2_6a[:, 2].reshape(b6a, 1)  # Trans meas force[mN] to [N]
    P6a = test2data6a['pressurePa']             # Intraluminal pressure [Pa]
    Lzz6a = test2data6a['lambda_z']             # Axial stretch
    ir6a = 1.0e-3*test2data6a['ir_mm']          # Inner radius [mm] to [m]
    or6a = 1.0e-3*test2data6a['or_mm']          # Outer radius [mm] to [m]
    
    h6a = or6a - ir6a                           # Thickness [m]
    Lqq6a = (ir6a + h6a/2.0)/(IR6 + H6/2.0)     # Circumferential stretch
    Lrr6a = 1.0/(Lzz6a*Lqq6a)                   # Radial stretch
    Tqq6a = (P6a*ir6a)/h6a                      # Circumferential Cauchy stress [Pa]
    Trr6a = (-P6a*ir6a)/(or6a + ir6a)           # Radial Cauchy stress [Pa]
    fp6a = np.pi*(ir6a**2)*P6a                  # Force due to pressure [N]
    f6a = fp6a + ft6a                           # Total axial force [N]
    Tzz6a = f6a/(np.pi*(or6a**2 - ir6a**2))     # Axial Cauchy stress [Pa]
    # =========================================================================
    #   Read in and format provided Matlab structure array
    data7 = sio.loadmat(file7, appendmat=True)
    for key, value in list(data7.items()):
        data7[key] = get_matlab_value(value)
    passive7 = data7['passive']
#    biaxial_Fl7 = passive7['biaxial_Fl']
    biaxial_Pd7 = passive7['biaxial_Pd']
    unloaded7 = passive7['unloaded']
    # =========================================================================
    #   Read in unloaded geometry
    OR7 = 1.0e-3*unloaded7['OR_mm']
    IR7 = 1.0e-3*unloaded7['IR_mm']
    H7 = OR7- IR7
    # =========================================================================
    #   Read in and format data for the in vivo test
    test2_7 = biaxial_Pd7['test2']
    test2data7 = biaxial_Pd7['test2data']
    OD_7 = test2_7[:, 1]                      # Outer diamter [um]
    b7 = test2_7[:, 2].size
    ft7 = 1.0e-3*test2_7[:, 2].reshape(b7, 1)  # Transducer measured force[mN] to [N]
    P7 = test2data7['pressurePa']            # Intraluminal pressure [Pa]
    Lzz7 = test2data7['lambda_z']            # Axial stretch
    ir7 = 1.0e-3*test2data7['ir_mm']         # Inner radius [mm] to [m]
    or7 = 1.0e-3*test2data7['or_mm']         # Outer radius [mm] to [m]
    
    h7 = or7 - ir7                       # Thickness [m]
    Lqq7 = (ir7 + h7/2.0)/(IR7 + H7/2.0)    # Circumferential stretch
    Lrr7 = 1.0/(Lzz7*Lqq7)                # Radial stretch
    Tqq7 = (P7*ir7)/h7                  # Circumferential Cauchy stress [Pa]
    Trr7 = (-P7*ir7)/(or7 + ir7)        # Radial Cauchy stress [Pa]
    fp7 = np.pi*(ir7**2)*P7              # Force due to pressure [N]
    f7 = fp7 + ft7                       # Total axial force [N]
    Tzz7 = f7/(np.pi*(or7**2 - ir7**2))  # Axial Cauchy stress [Pa]
    # =========================================================================
    #   Read in and format data for the 4% ABOVE test
    test2_7a = biaxial_Pd7['test3']
    test2data7a = biaxial_Pd7['test3data']
    OD_7a = test2_7a[:, 1]                      # Outer diamter [um]
    b7a = test2_7a[:, 2].size
    ft7a = 1.0e-3*test2_7a[:, 2].reshape(b7a, 1)  # Transducer measured force[mN] to [N]
    P7a = test2data7a['pressurePa']            # Intraluminal pressure [Pa]
    Lzz7a = test2data7a['lambda_z']            # Axial stretch
    ir7a = 1.0e-3*test2data7a['ir_mm']         # Inner radius [mm] to [m]
    or7a = 1.0e-3*test2data7a['or_mm']         # Outer radius [mm] to [m]
    
    h7a = or7a - ir7a                       # Thickness [m]
    Lqq7a = (ir7a + h7a/2.0)/(IR7 + H7/2.0)    # Circumferential stretch
    Lrr7a = 1.0/(Lzz7a*Lqq7a)                # Radial stretch
    Tqq7a = (P7a*ir7a)/h7a                  # Circumferential Cauchy stress [Pa]
    Trr7a = (-P7a*ir7a)/(or7a + ir7a)        # Radial Cauchy stress [Pa]
    fp7a = np.pi*(ir7a**2)*P7a              # Force due to pressure [N]
    f7a = fp7a + ft7a                       # Total axial force [N]
    Tzz7a = f7a/(np.pi*(or7a**2 - ir7a**2))  # Axial Cauchy stress [Pa]
    # =========================================================================
    #   Read in and format provided Matlab structure array
    data8 = sio.loadmat(file8, appendmat=True)
    for key, value in list(data8.items()):
        data8[key] = get_matlab_value(value)
    passive8 = data8['passive']
#    biaxial_Fl8 = passive8['biaxial_Fl']
    biaxial_Pd8 = passive8['biaxial_Pd']
    unloaded8 = passive8['unloaded']
    # =========================================================================
    #   Read in unloaded geometry
    OR8 = 1.0e-3*unloaded8['OR_mm']
    IR8 = 1.0e-3*unloaded8['IR_mm']
    H8 = OR8- IR8
    # =========================================================================
    #   Read in and format data for the in vivo test
    test2_8 = biaxial_Pd8['test2']
    test2data8 = biaxial_Pd8['test2data']
    OD_8 = test2_8[:, 1]                      # Outer diamter [um]
    b8 = test2_8[:, 2].size
    ft8 = 1.0e-3*test2_8[:, 2].reshape(b8, 1)  # Transducer measured force[mN] to [N]
    P8 = test2data8['pressurePa']            # Intraluminal pressure [Pa]
    Lzz8 = test2data8['lambda_z']            # Axial stretch
    ir8 = 1.0e-3*test2data8['ir_mm']         # Inner radius [mm] to [m]
    or8 = 1.0e-3*test2data8['or_mm']         # Outer radius [mm] to [m]
    
    h8 = or8 - ir8                       # Thickness [m]
    Lqq8 = (ir8 + h8/2.0)/(IR8 + H8/2.0)    # Circumferential stretch
    Lrr8 = 1.0/(Lzz8*Lqq8)                # Radial stretch
    Tqq8 = (P8*ir8)/h8                  # Circumferential Cauchy stress [Pa]
    Trr8 = (-P8*ir8)/(or8 + ir8)        # Radial Cauchy stress [Pa]
    fp8 = np.pi*(ir8**2)*P8              # Force due to pressure [N]
    f8 = fp8 + ft8                       # Total axial force [N]
    Tzz8 = f8/(np.pi*(or8**2 - ir8**2))  # Axial Cauchy stress [Pa]
    # =========================================================================
    #   Read in and format data for the 4% ABOVE test
    test2_8a = biaxial_Pd8['test3']
    test2data8a = biaxial_Pd8['test3data']
    OD_8a = test2_8a[:, 1]                      # Outer diamter [um]
    b8a = test2_8a[:, 2].size
    ft8a = 1.0e-3*test2_8a[:, 2].reshape(b8a, 1)  # Transducer measured force[mN] to [N]
    P8a = test2data8a['pressurePa']            # Intraluminal pressure [Pa]
    Lzz8a = test2data8a['lambda_z']            # Axial stretch
    ir8a = 1.0e-3*test2data8a['ir_mm']         # Inner radius [mm] to [m]
    or8a = 1.0e-3*test2data8a['or_mm']         # Outer radius [mm] to [m]
    
    h8a = or8a - ir8a                       # Thickness [m]
    Lqq8a = (ir8a + h8a/2.0)/(IR8 + H8/2.0)    # Circumferential stretch
    Lrr8a = 1.0/(Lzz8a*Lqq8a)                # Radial stretch
    Tqq8a = (P8a*ir8a)/h8a                  # Circumferential Cauchy stress [Pa]
    Trr8a = (-P8a*ir8a)/(or8a + ir8a)        # Radial Cauchy stress [Pa]
    fp8a = np.pi*(ir8a**2)*P8a              # Force due to pressure [N]
    f8a = fp8a + ft8a                       # Total axial force [N]
    Tzz8a = f8a/(np.pi*(or8a**2 - ir8a**2))  # Axial Cauchy stress [Pa]
    # =========================================================================
    #   Choose in vivo data only for fit
    OD_iv = np.mean([OD_1, OD_2, OD_3, OD_4, OD_5, OD_6, OD_7, OD_8], axis=0)
    OD_iv_std = np.std([OD_1, OD_2, OD_3, OD_4, OD_5, OD_6, OD_7, OD_8], axis=0)
    OD_iv_sem = sem([OD_1, OD_2, OD_3, OD_4, OD_5, OD_6, OD_7, OD_8], axis=0)
    
    or_iv = np.mean([or1, or2, or3, or4, or5, or6, or7, or8], axis=0)
    or_iv_std = np.std([or1, or2, or3, or4, or5, or6, or7, or8], axis=0)
    or_iv_sem = sem([or1, or2, or3, or4, or5, or6, or7, or8], axis=0)
    
    ir_iv = np.mean([ir1, ir2, ir3, ir4, ir5, ir6, ir7, ir8], axis=0)
    ir_iv_std = np.std([ir1, ir2, ir3, ir4, ir5, ir6, ir7, ir8], axis=0)
    ir_iv_sem = sem([ir1, ir2, ir3, ir4, ir5, ir6, ir7, ir8], axis=0)
    
    h_iv = np.mean([h1, h2, h3, h4, h5, h6, h7, h8], axis=0)
    h_iv_std = np.std([h1, h2, h3, h4, h5, h6, h7, h8], axis=0)
    h_iv_sem = sem([h1, h2, h3, h4, h5, h6, h7, h8], axis=0)
    
    Lzz_iv = np.mean([Lzz1, Lzz2, Lzz3, Lzz4, Lzz5, Lzz6, Lzz7, Lzz8], axis=0)
    Lzz_iv_std = np.std([Lzz1, Lzz2, Lzz3, Lzz4, Lzz5, Lzz6, Lzz7, Lzz8], axis=0)
    Lzz_iv_sem = sem([Lzz1, Lzz2, Lzz3, Lzz4, Lzz5, Lzz6, Lzz7, Lzz8], axis=0)
    
    Lqq_iv = np.mean([Lqq1, Lqq2, Lqq3, Lqq4, Lqq5, Lqq6, Lqq7, Lqq8], axis=0)
    Lqq_iv_std = np.std([Lqq1, Lqq2, Lqq3, Lqq4, Lqq5, Lqq6, Lqq7, Lqq8], axis=0)
    Lqq_iv_sem = sem([Lqq1, Lqq2, Lqq3, Lqq4, Lqq5, Lqq6, Lqq7, Lqq8], axis=0)
    
    Lrr_iv = np.mean([Lrr1, Lrr2, Lrr3, Lrr4, Lrr5, Lrr6, Lrr7, Lrr8], axis=0)
    Lrr_iv_std = np.std([Lrr1, Lrr2, Lrr3, Lrr4, Lrr5, Lrr6, Lrr7, Lrr8], axis=0)
    Lrr_iv_sem = sem([Lrr1, Lrr2, Lrr3, Lrr4, Lrr5, Lrr6, Lrr7, Lrr8], axis=0)
    
    P_iv = np.mean([P1, P2, P3, P4, P5, P6, P7, P8], axis=0)
    P_iv_std = np.std([P1, P2, P3, P4, P5, P6, P7, P8], axis=0)
    P_iv_sem = sem([P1, P2, P3, P4, P5, P6, P7, P8], axis=0)
    
    ft_iv = np.mean([ft1, ft2, ft3, ft4, ft5, ft6, ft7, ft8], axis=0)
    ft_iv_std = np.std([ft1, ft2, ft3, ft4, ft5, ft6, ft7, ft8], axis=0)
    ft_iv_sem = sem([ft1, ft2, ft3, ft4, ft5, ft6, ft7, ft8], axis=0)
    
    fp_iv = np.mean([fp1, fp2, fp3, fp4, fp5, fp6, fp7, fp8], axis=0)
    fp_iv_std = np.std([fp1, fp2, fp3, fp4, fp5, fp6, fp7, fp8], axis=0)
    fp_iv_sem = sem([fp1, fp2, fp3, fp4, fp5, fp6, fp7, fp8], axis=0)
    
    f_iv = np.mean([f1, f2, f3, f4, f5, f6, f7, f8], axis=0)
    f_iv_std = np.std([f1, f2, f3, f4, f5, f6, f7, f8], axis=0)
    f_iv_sem = sem([f1, f2, f3, f4, f5, f6, f7, f8], axis=0)
    
    Trr_iv = np.mean([Trr1, Trr2, Trr3, Trr4, Trr5, Trr6, Trr7, Trr8], axis=0)
    Trr_iv_std = np.std([Trr1, Trr2, Trr3, Trr4, Trr5, Trr6, Trr7, Trr8], axis=0)
    Trr_iv_sem = sem([Trr1, Trr2, Trr3, Trr4, Trr5, Trr6, Trr7, Trr8], axis=0)
    
    Tqq_iv = np.mean([Tqq1, Tqq2, Tqq3, Tqq4, Tqq5, Tqq6, Tqq7, Tqq8], axis=0)
    Tqq_iv_std = np.std([Tqq1, Tqq2, Tqq3, Tqq4, Tqq5, Tqq6, Tqq7, Tqq8], axis=0)
    Tqq_iv_sem = sem([Tqq1, Tqq2, Tqq3, Tqq4, Tqq5, Tqq6, Tqq7, Tqq8], axis=0)
    
    Tzz_iv = np.mean([Tzz1, Tzz2, Tzz3, Tzz4, Tzz5, Tzz6, Tzz7, Tzz8], axis=0)
    Tzz_iv_std = np.std([Tzz1, Tzz2, Tzz3, Tzz4, Tzz5, Tzz6, Tzz7, Tzz8], axis=0)
    Tzz_iv_sem = sem([Tzz1, Tzz2, Tzz3, Tzz4, Tzz5, Tzz6, Tzz7, Tzz8], axis=0)
    # =========================================================================
    #   4% ABOVE concatenated data
    OD_4a = np.mean([OD_1a, OD_2a, OD_3a, OD_4a, OD_5a, OD_6a, OD_7a, OD_8a], axis=0)
    OD_4a_std = np.std([OD_1a, OD_2a, OD_3a, OD_4a, OD_5a, OD_6a, OD_7a, OD_8a], axis=0)
    OD_4a_sem = sem([OD_1a, OD_2a, OD_3a, OD_4a, OD_5a, OD_6a, OD_7a, OD_8a], axis=0)
    
    or_4a = np.mean([or1a, or2a, or3a, or4a, or5a, or6a, or7a, or8a], axis=0)
    or_4a_std = np.std([or1a, or2a, or3a, or4a, or5a, or6a, or7a, or8a], axis=0)
    or_4a_sem = sem([or1a, or2a, or3a, or4a, or5a, or6a, or7a, or8a], axis=0)
    
    ir_4a = np.mean([ir1a, ir2a, ir3a, ir4a, ir5a, ir6a, ir7a, ir8a], axis=0)
    ir_4a_std = np.std([ir1a, ir2a, ir3a, ir4a, ir5a, ir6a, ir7a, ir8a], axis=0)
    ir_4a_sem = sem([ir1a, ir2a, ir3a, ir4a, ir5a, ir6a, ir7a, ir8a], axis=0)
    
    h_4a = np.mean([h1a, h2a, h3a, h4a, h5a, h6a, h7a, h8a], axis=0)
    h_4a_std = np.std([h1a, h2a, h3a, h4a, h5a, h6a, h7a, h8a], axis=0)
    h_4a_sem = sem([h1a, h2a, h3a, h4a, h5a, h6a, h7a, h8a], axis=0)
    
    Lzz_4a = np.mean([Lzz1a, Lzz2a, Lzz3a, Lzz4a, Lzz5a, Lzz6a, Lzz7a, Lzz8a], axis=0)
    Lzz_4a_std = np.std([Lzz1a, Lzz2a, Lzz3a, Lzz4a, Lzz5a, Lzz6a, Lzz7a, Lzz8a], axis=0)
    Lzz_4a_sem = sem([Lzz1a, Lzz2a, Lzz3a, Lzz4a, Lzz5a, Lzz6a, Lzz7a, Lzz8a], axis=0)
    
    Lqq_4a = np.mean([Lqq1a, Lqq2a, Lqq3a, Lqq4a, Lqq5a, Lqq6a, Lqq7a, Lqq8a], axis=0)
    Lqq_4a_std = np.std([Lqq1a, Lqq2a, Lqq3a, Lqq4a, Lqq5a, Lqq6a, Lqq7a, Lqq8a], axis=0)
    Lqq_4a_sem = sem([Lqq1a, Lqq2a, Lqq3a, Lqq4a, Lqq5a, Lqq6a, Lqq7a, Lqq8a], axis=0)
    
    Lrr_4a = np.mean([Lrr1a, Lrr2a, Lrr3a, Lrr4a, Lrr5a, Lrr6a, Lrr7a, Lrr8a], axis=0)
    Lrr_4a_std = np.std([Lrr1a, Lrr2a, Lrr3a, Lrr4a, Lrr5a, Lrr6a, Lrr7a, Lrr8a], axis=0)
    Lrr_4a_sem = sem([Lrr1a, Lrr2a, Lrr3a, Lrr4a, Lrr5a, Lrr6a, Lrr7a, Lrr8a], axis=0)
    
    P_4a = np.mean([P1a, P2a, P3a, P4a, P5a, P6a, P7a, P8a], axis=0)
    P_4a_std = np.std([P1a, P2a, P3a, P4a, P5a, P6a, P7a, P8a], axis=0)
    P_4a_sem = sem([P1a, P2a, P3a, P4a, P5a, P6a, P7a, P8a], axis=0)
    
    ft_4a = np.mean([ft1a, ft2a, ft3a, ft4a, ft5a, ft6a, ft7a, ft8a], axis=0)
    ft_4a_std = np.std([ft1a, ft2a, ft3a, ft4a, ft5a, ft6a, ft7a, ft8a], axis=0)
    ft_4a_sem = sem([ft1a, ft2a, ft3a, ft4a, ft5a, ft6a, ft7a, ft8a], axis=0)
    
    fp_4a = np.mean([fp1a, fp2a, fp3a, fp4a, fp5a, fp6a, fp7a, fp8a], axis=0)
    fp_4a_std = np.std([fp1a, fp2a, fp3a, fp4a, fp5a, fp6a, fp7a, fp8a], axis=0)
    fp_4a_sem = sem([fp1a, fp2a, fp3a, fp4a, fp5a, fp6a, fp7a, fp8a], axis=0)
    
    f_4a = np.mean([f1a, f2a, f3a, f4a, f5a, f6a, f7a, f8a], axis=0)
    f_4a_std = np.std([f1a, f2a, f3a, f4a, f5a, f6a, f7a, f8a], axis=0)
    f_4a_sem = sem([f1a, f2a, f3a, f4a, f5a, f6a, f7a, f8a], axis=0)
    
    Trr_4a = np.mean([Trr1a, Trr2a, Trr3a, Trr4a, Trr5a, Trr6a, Trr7a, Trr8a], axis=0)
    Trr_4a_std = np.std([Trr1a, Trr2a, Trr3a, Trr4a, Trr5a, Trr6a, Trr7a, Trr8a], axis=0)
    Trr_4a_sem = sem([Trr1a, Trr2a, Trr3a, Trr4a, Trr5a, Trr6a, Trr7a, Trr8a], axis=0)
    
    Tqq_4a = np.mean([Tqq1a, Tqq2a, Tqq3a, Tqq4a, Tqq5a, Tqq6a, Tqq7a, Tqq8a], axis=0)
    Tqq_4a_std = np.std([Tqq1a, Tqq2a, Tqq3a, Tqq4a, Tqq5a, Tqq6a, Tqq7a, Tqq8a], axis=0)
    Tqq_4a_sem = sem([Tqq1a, Tqq2a, Tqq3a, Tqq4a, Tqq5a, Tqq6a, Tqq7a, Tqq8a], axis=0)
    
    Tzz_4a = np.mean([Tzz1a, Tzz2a, Tzz3a, Tzz4a, Tzz5a, Tzz6a, Tzz7a, Tzz8a], axis=0)
    Tzz_4a_std = np.std([Tzz1a, Tzz2a, Tzz3a, Tzz4a, Tzz5a, Tzz6a, Tzz7a, Tzz8a], axis=0)
    Tzz_4a_sem = sem([Tzz1a, Tzz2a, Tzz3a, Tzz4a, Tzz5a, Tzz6a, Tzz7a, Tzz8a], axis=0)
    # =============================================================================
    #   Concatenation of in vivo with 4% above
    
    #Or = np.concatenate((or_iv, or_4a), axis=0)
    #Or.sort(axis=0)
    #Ir = np.concatenate((ir_iv, ir_4a), axis=0)
    #Ir.sort(axis=0)
    #h = np.concatenate((h_iv, h_4a), axis=0)
    #h.sort(axis=0)
    #Lzz = np.concatenate((Lzz_iv, Lzz_4a), axis=0)
    #Lzz.sort(axis=0)
    #Lqq = np.concatenate((Lqq_iv, Lqq_4a), axis=0)
    #Lqq.sort(axis=0)
    #Lrr = np.concatenate((Lrr_iv, Lrr_4a), axis=0)
    #Lrr.sort(axis=0)
    #P_exp = np.concatenate((P_iv, P_4a), axis=0)
    #P_exp.sort(axis=0)
    #fp_exp = np.concatenate((fp_iv, fp_4a), axis=0)
    #fp_exp.sort(axis=0)
    #ft_exp = np.concatenate((ft_iv, ft_4a), axis=0)
    #ft_exp.sort(axis=0)
    #f_exp = np.concatenate((f_iv, f_4a), axis=0)
    #f_exp.sort(axis=0)
    #Tqq_exp = np.concatenate((Tqq_iv, Tqq_4a), axis=0)
    #Tqq_exp.sort(axis=0)
    #Tzz_exp = np.concatenate((Tzz_iv, Tzz_4a), axis=0)
    #Tzz_exp.sort(axis=0)
    # =========================================================================
    #   Choose in vivo data only for fit
    Or = or_iv
    Ir = ir_iv
    h = h_iv
    Lzz = Lzz_iv
    Lqq = Lqq_iv
    Lrr = Lrr_iv
    P_exp = P_iv/1000 # converted to kPa to avoid overflow in np.exp 01/20/18
    ft_exp = ft_iv/1000 # converted to kN to avoid overflow in np.exp 01/20/18
    fp_exp = fp_iv
    f_exp = f_iv
    Tqq_exp = Tqq_iv
    Tzz_exp = Tzz_iv
    # =========================================================================
    #   Choose 4% above data only for fit
    #Or = or_4a
    #Ir = ir_4a
    #h = h_4a
    #Lzz = Lzz_4a
    #Lqq = Lqq_4a
    #Lrr = Lrr_4a
    #P_exp = P_4a/1000 # converted to kPa to avoid overflow in np.exp 01/20/18
    #fp_exp = fp_4a/1000 # converted to kN to avoid overflow in np.exp 01/20/18
    #ft_exp = ft_4a
    #f_exp = f_4a
    #Tqq_exp = Tqq_4a
    #Tzz_exp = Tzz_4a
    # =========================================================================
    #   fit to concatenated/or not data
    Crr = Lrr**2
    Cqq = Lqq**2
    Czz = Lzz**2
    
    ''' Save results to MATLAB .mat file '''
    # =========================================================================
    mdict = {'h_iv': np.reshape(h_iv,(12,1)), 'or_iv': np.reshape(or_iv,(12,1)),
             'ir_iv': np.reshape(ir_iv,(12,1)), 'OD_iv': np.reshape(OD_iv,(12,1)),
             'OD_iv_sem': np.reshape(OD_iv_sem,(12,1)), 'Lzz_iv': Lzz_iv, 
             'Lzz_iv_sem': Lzz_iv_sem, 'Lqq_iv': Lqq_iv, 'Lqq_iv_sem': Lqq_iv_sem,
             'Lrr_iv': Lrr_iv, 'Lrr_iv_sem': Lrr_iv_sem, 'P_iv': P_iv, 'P_iv_sem': P_iv_sem,
             'ft_iv': ft_iv, 'ft_iv_sem': ft_iv_sem, 'f_iv': f_iv, 'f_iv_sem': f_iv_sem, 
             'Trr_iv': Trr_iv, 'Trr_iv_sem': Trr_iv_sem, 'Tqq_iv': Tqq_iv, 
             'Tqq_iv_sem': Tqq_iv_sem, 'Tzz_iv': Tzz_iv, 'Tzz_iv_sem': Tzz_iv_sem,
             'OD_4a': np.reshape(OD_4a,(12,1)), 'OD_4a_sem': np.reshape(OD_4a_sem,(12,1)), 'Lzz_4a': Lzz_4a, 
             'Lzz_4a_sem': Lzz_4a_sem, 'Lqq_4a': Lqq_4a, 'Lqq_4a_sem': Lqq_4a_sem,
             'Lrr_4a': Lrr_4a, 'Lrr_4a_sem': Lrr_4a_sem, 'P_4a': P_4a, 'P_4a_sem': P_4a_sem,
             'ft_4a': ft_4a, 'ft_4a_sem': ft_4a_sem, 'f_4a': f_4a, 'f_4a_sem': f_4a_sem, 
             'Trr_4a': Trr_4a, 'Trr_4a_sem': Trr_4a_sem, 'Tqq_4a': Tqq_4a, 
             'Tqq_4a_sem': Tqq_4a_sem, 'Tzz_4a': Tzz_4a, 'Tzz_4a_sem': Tzz_4a_sem}
    
    sio.savemat("%s\%s\%s\%s" % (path, filename, model, 'exp_data'), mdict=mdict,
                appendmat=True)
    # =========================================================================

    return Or, Ir, h, Crr, Cqq, Czz, P_exp, ft_exp, Lrr_iv, Lqq_iv, Lzz_iv, ir_iv, or_iv, h_iv, ft_iv, f_iv, OD_iv, P_iv, Trr_iv, Tqq_iv, Tzz_iv, Lrr_4a, Lqq_4a, Lzz_4a, ir_4a, or_4a, h_4a, ft_4a, f_4a, OD_4a, P_4a, Trr_4a, Tqq_4a, Tzz_4a

