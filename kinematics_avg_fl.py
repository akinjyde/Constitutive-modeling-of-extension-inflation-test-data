# -*- coding: utf-8 -*-
''' This module computes the average kinematic quantities for the F-L test data '''

# Import needed modules
import scipy.io as sio
#from spyder_kernels.utils import iofuncs
#from spyder.utils.iofuncs import get_matlab_value
import scipy.io as sio
import numpy as np
import os
import matplotlib.ticker as tk
import matplotlib.pyplot as plt
from scipy.stats import sem


filename1 = 'SS091916KR-VAGL'
filename2 = 'SS091216KR-VAGL'
filename3 = 'SS082916KR-VAGL'
filename4 = 'SS081616KR-VAGL'
filename5 = 'SS121317KR-VAGL'
filename6 = 'SS011218KR-VAGL'
filename7 = 'SS012318KR-VAGL'
filename8 = 'SS012418KR-VAGL'

# Choose test: 8 mmHg (2) or 17 mmHg (3)

mm_a = 'test2'
mm_b = 'test2data' 

Pa2Hg = 0.00750062
Pa2kPa = 1e-3
kPa2Hg = 7.500617
kN2mN = 1e6
N2mN = 1e3

datafilepath = 'D:\Inflation_Extension_Modeling\Vaginal digestion\Data'
file1 = "%s\%s" % (datafilepath, filename1)
file2 = "%s\%s" % (datafilepath, filename2)
file3 = "%s\%s" % (datafilepath, filename3)
file4 = "%s\%s" % (datafilepath, filename4)
file5 = "%s\%s" % (datafilepath, filename5)
file6 = "%s\%s" % (datafilepath, filename6)
file7 = "%s\%s" % (datafilepath, filename7)
file8 = "%s\%s" % (datafilepath, filename8)

#def main(file1, file2, file3, file4, file5, file6, file7, file8, path, filename, model):
    
# Read in and format provided Matlab structure array
data1 = sio.loadmat(file1, appendmat=True)
for key, value in list(data1.items()):
    data1[key] = get_matlab_value(value)
passive1 = data1['passive']
biaxial_Fl1 = passive1['biaxial_Fl']
#    biaxial_Pd1 = passive1['biaxial_Pd']
unloaded1 = passive1['unloaded']
# =========================================================================
#   Read in unloaded geometry
OR1 = 1e-3*unloaded1['OR_mm']             # Outer (unloaded) radius [m]
IR1 = 1e-3*unloaded1['IR_mm']             # Inner (unloaded) radius [m]
H1 = OR1 - IR1                            # Thickness (unloaded) [m]
# =========================================================================
#   Read in and format data for the 17 mmHg test
test1 = biaxial_Fl1[mm_a]
P1_mmhg = test1[:,0]                      # Intraluminal pressure [mmHg]
OD_1 = test1[:, 1]                        # Outer diamter [um]
b1 = test1[:, 2].size
ft1 = 1e-3*test1[:, 2].reshape(b1, 1)     # Trans meas force[mN] to [N]

test1data = biaxial_Fl1[mm_b]
OD1_mm = test1data['od_mm']               # Outer diamter [mm]
P1 = test1data['pressurePa']              # Intraluminal pressure [Pa]
Lzz1 = test1data['lambda_z']              # Axial stretch
ir1 = 1e-3*test1data['ir_mm']             # Inner radius [mm] to [m]
or1 = 1e-3*test1data['or_mm']             # Outer radius [mm] to [m]

h1 = or1 - ir1                              # Thickness [m]
Lqq1 = (ir1 + h1/2.0)/(IR1 + H1/2.0)        # Circumferential stretch
Lrr1 = 1.0/(Lzz1*Lqq1)                      # Radial stretch
Tqq1 = (P1*ir1)/h1                          # Circ Cauchy stress [Pa]
Trr1 = (-P1*ir1)/(or1 + ir1)                # Radial Cauchy stress [Pa]
fp1 = np.pi*(ir1**2)*P1                     # Force due to pressure [N]
f1 = fp1 + ft1                              # Total axial force [N]
Tzz1 = f1/(np.pi*(or1**2 - ir1**2))         # Axial Cauchy stress [Pa]
    # =========================================================================
#   Read in and format provided Matlab structure array
data2 = sio.loadmat(file2, appendmat=True)
for key, value in list(data2.items()):
    data2[key] = get_matlab_value(value)
passive2 = data2['passive']
biaxial_Fl2 = passive2['biaxial_Fl']
#    biaxial_Pd2 = passive2['biaxial_Pd']
unloaded2 = passive2['unloaded']
# =========================================================================
#   Read in unloaded geometry
OR2 = 1.0e-3*unloaded2['OR_mm']
IR2 = 1.0e-3*unloaded2['IR_mm']
H2 = OR2 - IR2
# =========================================================================
#   Read in and format data for the 8 mmHg test
test2 = biaxial_Fl2[mm_a]
P2_mmhg = test2[:,0]                     # Intraluminal pressure [mmHg]
OD_2 = test2[:, 1]                        # Outer diamter [um]
b2 = test2[:, 2].size
ft2 = 1e-3*test2[:, 2].reshape(b2, 1)   # Trans meas force[mN] to [N]

test2data = biaxial_Fl2[mm_b]
OD2_mm = test2data['od_mm']                # Outer diamter [mm]
P2 = test2data['pressurePa']            # Intraluminal pressure [Pa]
Lzz2 = test2data['lambda_z']               # Axial stretch
ir2 = 1e-3*test2data['ir_mm']            # Inner radius [mm] to [m]
or2 = 1e-3*test2data['or_mm']            # Outer radius [mm] to [m]

h2 = or2 - ir2                              # Thickness [m]
Lqq2 = (ir2 + h2/2.0)/(IR2 + H2/2.0)        # Circumferential stretch
Lrr2 = 1.0/(Lzz2*Lqq2)                      # Radial stretch
Tqq2 = (P2*ir2)/h2                          # Circ Cauchy stress [Pa]
Trr2 = (-P2*ir2)/(or2 + ir2)                # Radial Cauchy stress [Pa]
fp2 = np.pi*(ir2**2)*P2                     # Force due to pressure [N]
f2 = fp2 + ft2                              # Total axial force [N]
Tzz2 = f2/(np.pi*(or2**2 - ir2**2))         # Axial Cauchy stress [Pa]
# =========================================================================
#   Read in and format provided Matlab structure array
data3 = sio.loadmat(file3, appendmat=True)
for key, value in list(data3.items()):
    data3[key] = get_matlab_value(value)
passive3= data3['passive']
biaxial_Fl3= passive3['biaxial_Fl']
#    biaxial_Pd3= passive3['biaxial_Pd']
unloaded3= passive3['unloaded']
# =========================================================================
#   Read in unloaded geometry
OR3= 1e-3*unloaded3['OR_mm']
IR3= 1e-3*unloaded3['IR_mm']
H3= OR3- IR3
# =========================================================================
#   Read in and format data for the 8 mmHg test
test3 = biaxial_Fl3[mm_a]
P3_mmhg = test3[:,0]                     # Intraluminal pressure [mmHg]
OD_3 = test3[:, 1]                        # Outer diamter [um]
b3 = test3[:, 2].size
ft3 = 1e-3*test3[:, 2].reshape(b3, 1)   # Trans meas force[mN] to [N]

test3data = biaxial_Fl3[mm_b]
OD3_mm = test3data['od_mm']                # Outer diamter [mm]
P3 = test3data['pressurePa']            # Intraluminal pressure [Pa]
Lzz3 = test3data['lambda_z']               # Axial stretch
ir3 = 1e-3*test3data['ir_mm']            # Inner radius [mm] to [m]
or3 = 1e-3*test3data['or_mm']            # Outer radius [mm] to [m]

h3 = or3 - ir3                              # Thickness [m]
Lqq3 = (ir3 + h3/2.0)/(IR3 + H3/2.0)        # Circumferential stretch
Lrr3 = 1.0/(Lzz3*Lqq3)                      # Radial stretch
Tqq3 = (P3*ir3)/h3                          # Circ Cauchy stress [Pa]
Trr3 = (-P3*ir3)/(or3 + ir3)                # Radial Cauchy stress [Pa]
fp3 = np.pi*(ir3**2)*P3                     # Force due to pressure [N]
f3 = fp3 + ft3                              # Total axial force [N]
Tzz3 = f3/(np.pi*(or3**2 - ir3**2))         # Axial Cauchy stress [Pa]
# =========================================================================
#   Read in and format provided Matlab structure array
data4 = sio.loadmat(file4, appendmat=True)
for key, value in list(data4.items()):
    data4[key] = get_matlab_value(value)
passive4 = data4['passive']
biaxial_Fl4 = passive4['biaxial_Fl']
#    biaxial_Pd4 = passive4['biaxial_Pd']
unloaded4 = passive4['unloaded']
# =========================================================================
#   Read in unloaded geometry
OR4 = 1.0e-3*unloaded4['OR_mm']
IR4 = 1.0e-3*unloaded4['IR_mm']
H4 = OR4- IR4
# =========================================================================
#   Read in and format data for the 8 mmHg test
test4 = biaxial_Fl4[mm_a]
P4_mmhg = test4[:,0]                     # Intraluminal pressure [mmHg]
OD_4 = test4[:, 1]                        # Outer diamter [um]
b4 = test4[:, 2].size
ft4 = 1e-3*test4[:, 2].reshape(b4, 1)   # Trans meas force[mN] to [N]

test4data = biaxial_Fl4[mm_b]
OD4_mm = test4data['od_mm']                # Outer diamter [mm]
P4 = test4data['pressurePa']            # Intraluminal pressure [Pa]
Lzz4 = test4data['lambda_z']               # Axial stretch
ir4 = 1e-3*test4data['ir_mm']            # Inner radius [mm] to [m]
or4 = 1e-3*test4data['or_mm']            # Outer radius [mm] to [m]

h4 = or4 - ir4                              # Thickness [m]
Lqq4 = (ir4 + h4/2.0)/(IR4 + H4/2.0)        # Circumferential stretch
Lrr4 = 1.0/(Lzz4*Lqq4)                      # Radial stretch
Tqq4 = (P4*ir4)/h4                          # Circ Cauchy stress [Pa]
Trr4 = (-P4*ir4)/(or4 + ir4)                # Radial Cauchy stress [Pa]
fp4 = np.pi*(ir4**2)*P4                     # Force due to pressure [N]
f4 = fp4 + ft4                              # Total axial force [N]
Tzz4 = f4/(np.pi*(or4**2 - ir4**2))         # Axial Cauchy stress [Pa]
# =========================================================================
#   Read in and format provided Matlab structure array
data5 = sio.loadmat(file5, appendmat=True)
for key, value in list(data5.items()):
    data5[key] = get_matlab_value(value)
passive5 = data5['passive']
biaxial_Fl5 = passive5['biaxial_Fl']
#    biaxial_Pd5 = passive5['biaxial_Pd']
unloaded5 = passive5['unloaded']
# =========================================================================
#   Read in unloaded geometry
OR5 = 1.0e-3*unloaded5['OR_mm']
IR5 = 1.0e-3*unloaded5['IR_mm']
H5 = OR5- IR5
# =========================================================================
#   Read in and format data for the 8 mmHg test
test5 = biaxial_Fl5[mm_a]
P5_mmhg = test5[:,0]                     # Intraluminal pressure [mmHg]
OD_5 = test5[:, 1]                        # Outer diamter [um]
b5 = test5[:, 2].size
ft5 = 1e-3*test5[:, 2].reshape(b5, 1)   # Trans meas force[mN] to [N]

test5data = biaxial_Fl5[mm_b]
OD5_mm = test5data['od_mm']                # Outer diamter [mm]
P5 = test5data['pressurePa']            # Intraluminal pressure [Pa]
Lzz5 = test5data['lambda_z']               # Axial stretch
ir5 = 1e-3*test5data['ir_mm']            # Inner radius [mm] to [m]
or5 = 1e-3*test5data['or_mm']            # Outer radius [mm] to [m]

h5 = or5 - ir5                              # Thickness [m]
Lqq5 = (ir5 + h5/2.0)/(IR5 + H5/2.0)        # Circumferential stretch
Lrr5 = 1.0/(Lzz5*Lqq5)                      # Radial stretch
Tqq5 = (P5*ir5)/h5                          # Circ Cauchy stress [Pa]
Trr5 = (-P5*ir5)/(or5 + ir5)                # Radial Cauchy stress [Pa]
fp5 = np.pi*(ir5**2)*P5                     # Force due to pressure [N]
f5 = fp5 + ft5                              # Total axial force [N]
Tzz5 = f5/(np.pi*(or5**2 - ir5**2))         # Axial Cauchy stress [Pa]
# =========================================================================
#   Read in and format provided Matlab structure array
data6 = sio.loadmat(file6, appendmat=True)
for key, value in list(data6.items()):
    data6[key] = get_matlab_value(value)
passive6 = data6['passive']
biaxial_Fl6 = passive6['biaxial_Fl']
#biaxial_Pd6 = passive6['biaxial_Pd']
unloaded6 = passive6['unloaded']
# =========================================================================
#   Read in unloaded geometry
OR6 = 1e-3*unloaded6['OR_mm']
IR6 = 1e-3*unloaded6['IR_mm']
H6 = OR6- IR6
# =========================================================================
#   Read in and format data for the 8 mmHg test
test6 = biaxial_Fl6[mm_a]
P6_mmhg = test6[:,0]                     # Intraluminal pressure [mmHg]
OD_6 = test6[:, 1]                        # Outer diamter [um]
b6 = test6[:, 2].size
ft6 = 1e-3*test6[:, 2].reshape(b6, 1)   # Trans meas force[mN] to [N]

test6data = biaxial_Fl6[mm_b]
OD5_mm = test6data['od_mm']                # Outer diamter [mm]
P6 = test6data['pressurePa']            # Intraluminal pressure [Pa]
Lzz6 = test6data['lambda_z']               # Axial stretch
ir6 = 1e-3*test6data['ir_mm']            # Inner radius [mm] to [m]
or6 = 1e-3*test6data['or_mm']            # Outer radius [mm] to [m]

h6 = or6 - ir6                              # Thickness [m]
Lqq6 = (ir6 + h6/2.0)/(IR6 + H6/2.0)        # Circumferential stretch
Lrr6 = 1.0/(Lzz6*Lqq6)                      # Radial stretch
Tqq6 = (P6*ir6)/h6                          # Circ Cauchy stress [Pa]
Trr6 = (-P6*ir6)/(or6 + ir6)                # Radial Cauchy stress [Pa]
fp6 = np.pi*(ir6**2)*P6                     # Force due to pressure [N]
f6 = fp6 + ft6                              # Total axial force [N]
Tzz6 = f6/(np.pi*(or6**2 - ir6**2))         # Axial Cauchy stress [Pa]
# =========================================================================
#   Read in and format provided Matlab structure array
data7 = sio.loadmat(file7, appendmat=True)
for key, value in list(data7.items()):
    data7[key] = get_matlab_value(value)
passive7 = data7['passive']
biaxial_Fl7 = passive7['biaxial_Fl']
#    biaxial_Pd7 = passive7['biaxial_Pd']
unloaded7 = passive7['unloaded']
# =========================================================================
#   Read in unloaded geometry
OR7 = 1.0e-3*unloaded7['OR_mm']
IR7 = 1.0e-3*unloaded7['IR_mm']
H7 = OR7- IR7
# =========================================================================
#   Read in and format data for the 8 mmHg test
test7 = biaxial_Fl7[mm_a]
P7_mmhg = test7[:,0]                     # Intraluminal pressure [mmHg]
OD_7 = test7[:, 1]                        # Outer diamter [um]
b7 = test7[:, 2].size
ft7 = 1.0e-3*test7[:, 2].reshape(b7, 1)   # Trans meas force[mN] to [N]

test7data = biaxial_Fl7[mm_b]
OD7_mm = test7data['od_mm']                # Outer diamter [mm]
P7 = test7data['pressurePa']            # Intraluminal pressure [Pa]
Lzz7 = test7data['lambda_z']               # Axial stretch
ir7 = 1.0e-3*test7data['ir_mm']            # Inner radius [mm] to [m]
or7 = 1.0e-3*test7data['or_mm']            # Outer radius [mm] to [m]

h7 = or7 - ir7                              # Thickness [m]
Lqq7 = (ir7 + h7/2.0)/(IR7 + H7/2.0)        # Circumferential stretch
Lrr7 = 1.0/(Lzz7*Lqq7)                      # Radial stretch
Tqq7 = (P7*ir7)/h7                          # Circ Cauchy stress [Pa]
Trr7 = (-P7*ir7)/(or7 + ir7)                # Radial Cauchy stress [Pa]
fp7 = np.pi*(ir7**2)*P7                     # Force due to pressure [N]
f7 = fp7 + ft7                              # Total axial force [N]
Tzz7 = f7/(np.pi*(or7**2 - ir7**2))         # Axial Cauchy stress [Pa]
# =========================================================================
#   Read in and format provided Matlab structure array
data8 = sio.loadmat(file8, appendmat=True)
for key, value in list(data8.items()):
    data8[key] = get_matlab_value(value)
passive8 = data8['passive']
biaxial_Fl8 = passive8['biaxial_Fl']
#    biaxial_Pd8 = passive8['biaxial_Pd']
unloaded8 = passive8['unloaded']
# =========================================================================
#   Read in unloaded geometry
OR8 = 1.0e-3*unloaded8['OR_mm']
IR8 = 1.0e-3*unloaded8['IR_mm']
H8 = OR8- IR8
# =========================================================================
#   Read in and format data for the 8 mmHg test
test8 = biaxial_Fl8[mm_a]
P8_mmhg = test8[:,0]                     # Intraluminal pressure [mmHg]
OD_8 = test8[:, 1]                        # Outer diamter [um]
b8 = test8[:, 2].size
ft8 = 1.0e-3*test8[:, 2].reshape(b8, 1)   # Trans meas force[mN] to [N]

test8data = biaxial_Fl8[mm_b]
OD8_mm = test8data['od_mm']                # Outer diamter [mm]
P8 = test8data['pressurePa']            # Intraluminal pressure [Pa]
Lzz8 = test8data['lambda_z']               # Axial stretch
ir8 = 1.0e-3*test8data['ir_mm']            # Inner radius [mm] to [m]
or8 = 1.0e-3*test8data['or_mm']            # Outer radius [mm] to [m]

h8 = or8 - ir8                              # Thickness [m]
Lqq8 = (ir8 + h8/2.0)/(IR8 + H8/2.0)        # Circumferential stretch
Lrr8 = 1.0/(Lzz8*Lqq8)                      # Radial stretch
Tqq8 = (P8*ir8)/h8                          # Circ Cauchy stress [Pa]
Trr8 = (-P8*ir8)/(or8 + ir8)                # Radial Cauchy stress [Pa]
fp8 = np.pi*(ir8**2)*P8                     # Force due to pressure [N]
f8 = fp8 + ft8                              # Total axial force [N]
Tzz8 = f8/(np.pi*(or8**2 - ir8**2))         # Axial Cauchy stress [Pa]
# =========================================================================
#   Choose in vivo data only for fit
OD = np.mean([OD_1, OD_2, OD_3, OD_4, OD_5, OD_6, OD_7, OD_8], axis=0)
OD_std = np.std([OD_1, OD_2, OD_3, OD_4, OD_5, OD_6, OD_7, OD_8], axis=0)
OD_sem = sem([OD_1, OD_2, OD_3, OD_4, OD_5, OD_6, OD_7, OD_8], axis=0)

Or = np.mean([or1, or2, or3, or4, or5, or6, or7, or8], axis=0)
Or_std = np.std([or1, or2, or3, or4, or5, or6, or7, or8], axis=0)
Or_sem = sem([or1, or2, or3, or4, or5, or6, or7, or8], axis=0)

Ir = np.mean([ir1, ir2, ir3, ir4, ir5, ir6, ir7, ir8], axis=0)
Ir_std = np.std([ir1, ir2, ir3, ir4, ir5, ir6, ir7, ir8], axis=0)
Ir_sem = sem([ir1, ir2, ir3, ir4, ir5, ir6, ir7, ir8], axis=0)

h = np.mean([h1, h2, h3, h4, h5, h6, h7, h8], axis=0)
h_std = np.std([h1, h2, h3, h4, h5, h6, h7, h8], axis=0)
h_sem = sem([h1, h2, h3, h4, h5, h6, h7, h8], axis=0)

Lzz = np.mean([Lzz1, Lzz2, Lzz3, Lzz4, Lzz5, Lzz6, Lzz7, Lzz8], axis=0)
Lzz_std = np.std([Lzz1, Lzz2, Lzz3, Lzz4, Lzz5, Lzz6, Lzz7, Lzz8], axis=0)
Lzz_sem = sem([Lzz1, Lzz2, Lzz3, Lzz4, Lzz5, Lzz6, Lzz7, Lzz8], axis=0)

Lqq = np.mean([Lqq1, Lqq2, Lqq3, Lqq4, Lqq5, Lqq6, Lqq7, Lqq8], axis=0)
Lqq_std = np.std([Lqq1, Lqq2, Lqq3, Lqq4, Lqq5, Lqq6, Lqq7, Lqq8], axis=0)
Lqq_sem = sem([Lqq1, Lqq2, Lqq3, Lqq4, Lqq5, Lqq6, Lqq7, Lqq8], axis=0)

Lrr = np.mean([Lrr1, Lrr2, Lrr3, Lrr4, Lrr5, Lrr6, Lrr7, Lrr8], axis=0)
Lrr_std = np.std([Lrr1, Lrr2, Lrr3, Lrr4, Lrr5, Lrr6, Lrr7, Lrr8], axis=0)
Lrr_sem = sem([Lrr1, Lrr2, Lrr3, Lrr4, Lrr5, Lrr6, Lrr7, Lrr8], axis=0)

P = np.mean([P1, P2, P3, P4, P5, P6, P7, P8], axis=0)
P_std = np.std([P1, P2, P3, P4, P5, P6, P7, P8], axis=0)
P_sem = sem([P1, P2, P3, P4, P5, P6, P7, P8], axis=0)

ft = np.mean([ft1, ft2, ft3, ft4, ft5, ft6, ft7, ft8], axis=0)
ft_std = np.std([ft1, ft2, ft3, ft4, ft5, ft6, ft7, ft8], axis=0)
ft_sem = sem([ft1, ft2, ft3, ft4, ft5, ft6, ft7, ft8], axis=0)

fp = np.mean([fp1, fp2, fp3, fp4, fp5, fp6, fp7, fp8], axis=0)
fp_std = np.std([fp1, fp2, fp3, fp4, fp5, fp6, fp7, fp8], axis=0)
fp_sem = sem([fp1, fp2, fp3, fp4, fp5, fp6, fp7, fp8], axis=0)

f = np.mean([f1, f2, f3, f4, f5, f6, f7, f8], axis=0)
f_std = np.std([f1, f2, f3, f4, f5, f6, f7, f8], axis=0)
f_sem = sem([f1, f2, f3, f4, f5, f6, f7, f8], axis=0)

Trr = np.mean([Trr1, Trr2, Trr3, Trr4, Trr5, Trr6, Trr7, Trr8], axis=0)
Trr_std = np.std([Trr1, Trr2, Trr3, Trr4, Trr5, Trr6, Trr7, Trr8], axis=0)
Trr_sem = sem([Trr1, Trr2, Trr3, Trr4, Trr5, Trr6, Trr7, Trr8], axis=0)

Tqq = np.mean([Tqq1, Tqq2, Tqq3, Tqq4, Tqq5, Tqq6, Tqq7, Tqq8], axis=0)
Tqq_std = np.std([Tqq1, Tqq2, Tqq3, Tqq4, Tqq5, Tqq6, Tqq7, Tqq8], axis=0)
Tqq_sem = sem([Tqq1, Tqq2, Tqq3, Tqq4, Tqq5, Tqq6, Tqq7, Tqq8], axis=0)

Tzz = np.mean([Tzz1, Tzz2, Tzz3, Tzz4, Tzz5, Tzz6, Tzz7, Tzz8], axis=0)
Tzz_std = np.std([Tzz1, Tzz2, Tzz3, Tzz4, Tzz5, Tzz6, Tzz7, Tzz8], axis=0)
Tzz_sem = sem([Tzz1, Tzz2, Tzz3, Tzz4, Tzz5, Tzz6, Tzz7, Tzz8], axis=0)

# =========================================================================
#   Choose in vivo data only for fit
P_exp = P/1000 # converted to kPa to avoid overflow in np.exp 01/20/18
ft_exp = ft/1000 # converted to kN to avoid overflow in np.exp 01/20/18
fp_exp = fp
f_exp = f
# =========================================================================
#   fit to concatenated/or not data
Crr = Lrr**2
Cqq = Lqq**2
Czz = Lzz**2

fs=12

#plt.figure()
#plt.plot(Lzz1, f1*1000, 'o-', lw=0.5, ms=3, mew=0.5, fillstyle='none', color='k', label='091916')
#plt.plot(Lzz2, f2*1000, 'd--', lw=0.5, ms=3, mew=0.5, fillstyle='none', color='r', label='091216')
#plt.plot(Lzz3, f3*1000, 's-', lw=0.5, ms=3, mew=0.5, fillstyle='none', color='b', label='082916')
#plt.plot(Lzz4, f4*1000, '^-', lw=0.5, ms=3, mew=0.5, fillstyle='none', color='y', label='081616')
#plt.plot(Lzz5, f5*1000, '+-', lw=0.5, ms=3, mew=0.5, fillstyle='none', color='c', label='121317')
#plt.plot(Lzz6, f6*1000, 'p-', lw=0.5, ms=3, mew=0.5, fillstyle='none', color='m', label='011218')
#plt.plot(Lzz7, f7*1000, 'v-', lw=0.5, ms=3, mew=0.5, fillstyle='none', color='g', label='012318')
#plt.plot(Lzz8, f8*1000, 'x-', lw=0.5, ms=3, mew=0.5, fillstyle='none', color='y', label='012418')
#plt.xlabel(r'$\lambda_z$', fontsize=fs)
#plt.ylabel("Axial Force (mN)", fontsize=fs)
#plt.legend(loc='upper left', fontsize=fs, frameon=False)
#plt.tick_params(axis='both', labelsize=fs, width=1, length=5)
#plt.xticks(rotation=45, fontsize=fs)
#plt.xlim(1.05, 1.28)
#plt.ylim(0, 80)
#ax = plt.gca()
#ax.xaxis.set_major_formatter(tk.FormatStrFormatter('%.2f'))
#ax.yaxis.set_major_locator(tk.MultipleLocator(10))
#ax.xaxis.set_major_locator(tk.MultipleLocator(0.02))
#ax.spines['right'].set_visible(False)
#ax.spines['top'].set_visible(False)
#plt.tight_layout()

path1 = r'C:\Users\Jyde\Desktop'

#plt.savefig("%s\%s" % (path1, 'F_L_plot1.png'), dpi=1200, format='png',
#            bbox_inches="tight")

plt.figure()
# Plot of experimental data (avg +/- sem) for 8mmHg F-L tests (n=8) 
l1=plt.errorbar(Lzz, Tzz*Pa2kPa, xerr=Lzz_sem, yerr=Tzz_sem*Pa2kPa, fmt='ok',
                fillstyle='none', label='control (@ 8 mmHg)', lw=1, elinewidth=0.5,
                capsize=1)

#l1=plt.errorbar(Lzz, ft_exp*1e6, xerr=Lzz_sem, yerr=ft_sem*1e3, fmt='ok',
#                fillstyle='none', label='control (@ 17 mmHg)', lw=1, elinewidth=0.5,
#                capsize=1)

#   Control
a = np.deg2rad(55.06557222)
mu = 2.807085817
c1 = 0.271978324
c2 = 63.01000321

#   Elastase
#a = np.deg2rad(55.9356105)
#mu = 2.262625512
#c1 = 0.18423886
#c2 = 102.8309591


I4 = Cqq*(np.sin(a))**2 + Czz*(np.cos(a))**2
W1 = 0.5*mu     # kPa
W4 = 2.0*c1*(I4 - 1.0)*np.exp(c2*(I4 - 1.0)**2)     # kPa
P_th = 2.0*(W1*(Cqq - Crr) + W4*Cqq*(np.sin(a))**2)*np.log(Or/Ir)  # kPa

f_int = 2.0*(W1*(2.0*Czz - Crr - Cqq) + W4*(2.0*Czz*(np.cos(a))**2 - Cqq*(np.sin(a))**2))   # kPa
ft_th = 0.5*np.pi*(Or**2 - Ir**2)*f_int                 # kN i.e. m^2 * kPa

fp_th = np.pi*(Ir**2)*P_exp                     # Force due to pressure [kN]
f_th = fp_th + ft_th                              # Total axial force [kN]
Tzz_th = f_th/(np.pi*(Or**2 - Ir**2))         # Axial Cauchy stress [kPa]



plt.plot(Lzz, Tzz_th)



plt.xlabel(r'$\lambda_z$', fontsize=fs)
plt.ylabel(r'$\sigma_{zz}$ (kPa)', fontsize=fs)
plt.legend(loc='upper left', fontsize=fs, frameon=False)
plt.tick_params(axis='both', labelsize=fs, width=1, length=5)
plt.xticks(rotation=0, fontsize=fs)

#plt.savefig("%s\%s" % (path1, 'F_L_ctl-P_exp_8.png'), dpi=1200, format='png',
#            bbox_inches="tight")
#plt.figure()
#plt.plot(Lzz, P_exp*1e3*Pa2Hg, label ='P_exp')
#plt.plot(Lzz, P_th*1e3*Pa2Hg, label = 'P_th')
#plt.legend(loc='upper left', fontsize=fs, frameon=False)
Tqq_th = (P_th*Ir)/h                   # Circ Cauchy stress [kPa]
plt.figure()
l1=plt.errorbar(Lzz, Tqq*Pa2kPa, xerr=Lzz_sem, yerr=Tqq_sem*Pa2kPa, fmt='ok',
                fillstyle='none', label='control (@ 8 mmHg)', lw=1, elinewidth=0.5,
                capsize=1)
plt.plot(Lzz, Tqq_th)

#plt.figure()
#plt.plot(Lzz, label='Lzz')
#plt.plot(Lqq, label='Lqq')
#plt.plot(Lrr, label='Lrr')
#plt.plot(Lzz*Lqq*Lrr, label='Lqq*Lzz*Lrr')
#plt.legend(loc='best', fontsize=fs, frameon=False)
# =============================================================================
"""Specify .mat filename """
# =============================================================================
filename = 'avg_ctl_load'
# =============================================================================
""" Specify model """
# =============================================================================
model = 'NH_2FF'
# =============================================================================
z = "%s\%s" % (filename, model)
if not os.path.exists(z):
    os.makedirs(z)
    
path = 'D:\Inflation_Extension_Modeling\Vaginal digestion'

''' Save experimental and theoretical Fl test data to MATLAB .mat file '''
# =============================================================================
#mdict = {'Tzz': Tzz, 'Tzz_sem': Tzz_sem, 'Lzz': Lzz, 'Lzz_sem': Lzz_sem,
#         'Tzz_th': Tzz_th}
#
#sio.savemat("%s\%s\%s\%s" % (path, filename, model, 'FL_data_ctl'), mdict=mdict,
#            appendmat=True)