# -*- coding: utf-8 -*-
"""
Created by Akinjide R. Akintunde 
January 2018 

Code fits inflation-extension experimental data for transversely-isotropic 
tissues to: NH_2FF i.e. neo-Hookean + 2 Fiber Family model
(Holzapfel et al., J Elast. 2000)
https://link.springer.com/article/10.1023/A:1010835316564

Note: This code makes use of spyder-kernels Input/Output module i.e. 
utils.iofuncs to read in data from .mat files through the function 
get_matlab_value(). Hence, the code is best run in Spyder IDE otherwise the 
dependency should be removed in kinematics.py
"""
# =============================================================================
# Import created modules
from kinematics import main as runKin
from solve4Lqq_iv import main as runSolve4Lqq_iv
from lin_stiff import main as runLinStiff
from sef_contour import main as runSEFContour
from corr_mat_sens_analysis import main as runCorrSens
from slope_intercept import main as runSlopeIntercept
from plot import main as runPlot
from exp_stiffness_compliance import main as runStiffComp
from fit_diff_evol import main as runDiffEvolFit


# Import needed standard Python modules
import numpy as np
import os
import xlsxwriter
# =============================================================================
""" Specify .mat filename """
# =============================================================================
filename = 'control1'

datafilepath = 'D:\murine-vagina\Data'
path = 'D:\murine-vagina'


file = "%s\%s" % (datafilepath, filename)
# =============================================================================
z = "%s" % (filename)
if not os.path.exists(z):
    os.makedirs(z)


(Or, Ir, h, Crr, Cqq, Czz, P_exp, ft_exp, Lrriv, Lqqiv, Lzziv, iriv, oriv, hiv,
 ftiv, fiv, OD_iv, P_iv, Trriv, Tqqiv, Tzziv, Lrr4a, Lqq4a, Lzz4a, ir4a, or4a,
 h4a, ft4a, f4a, OD_4a, P_4a, Trr4a, Tqq4a, Tzz4a) = runKin(file)

# Specify parameter bounds 
#              ce     c1     c2     alpha
#              kPa    kPa    -      rad
lb = np.array([0.0,   0.0,   0.0,   0.0])
ub = np.array([1.0e2, 1.0e2, 1.0e3, (np.pi/2.0)])    

bounds = (np.array([(lb[0], ub[0]), (lb[1], ub[1]),
                    (lb[2], ub[2]), (lb[3], ub[3])]))
    
# Perform optimization of obj. fun. using differential evolution algorithm
resDE, gof = runDiffEvolFit(bounds, Crr, Cqq, Czz, Ir, Or, P_exp, ft_exp,
                            filename, path)
x = resDE.x

# Solve for in vivo circumferential stretch
Lqq_iv, Lzz_iv, ir_iv, or_iv = runSolve4Lqq_iv(x, P_iv, Lqqiv, Lzziv, iriv,
                                               oriv, filename, path)

# Compute linearized stiffness
Crrrr_iv, Cqqqq_iv, Czzzz_iv = runLinStiff(x, Lqq_iv, Lzz_iv, filename, path)

# Compute the slope at the chosen in vivo circumferential stretch
dPdLq, PIntercept = runSlopeIntercept(x, Lqq_iv, Lzz_iv, ir_iv, or_iv, filename)

# Compute in vivo stored energy and make contour plot
Wiv = runSEFContour(x, Lqq_iv, Lzz_iv, filename, path)

# Plot and save figures
(Trrivth, Tqqivth, Tzzivth, Tqzivth, P_th_iv, ft_th_iv, f_th_iv, Trr4ath,
 Tqq4ath, Tzz4ath, P_th_4a, ft_th_4a, f_th_4a) = runPlot(x, Lrriv, Lqqiv, Lzziv,
                                             iriv, oriv, hiv, ftiv, fiv, OD_iv,
                                             P_iv, Trriv, Tqqiv, Tzziv, Lrr4a,
                                             Lqq4a, Lzz4a, ir4a, or4a, h4a,
                                             ft4a, f4a, OD_4a, P_4a, Trr4a,
                                             Tqq4a, Tzz4a, filename, path,
                                             dPdLq, PIntercept, Lqq_iv)

# Compute correlation matrix and sensitivity analysis
Corr_Mat_det = runCorrSens(x, Crr, Cqq, Czz, Or, Ir, h, Tqqivth, Tzzivth,
                           filename, path)

# Interpolate outer diameter values and compute compliance*
(OD_2_5, OD_5, OD_7_5, OD_10, OD_15, OD_20,
 C1, C2, C3, C4, C5) = runStiffComp(P_iv, OD_iv)

# =============================================================================
''' Save results to MS Excel file '''

# make dictionary
results = (['ce_kPa', resDE.x[0]], ['c1_kPa', resDE.x[1]],
           ['c2', resDE.x[2]], ['alpha_deg', np.degrees(resDE.x[3])],
           ['obj_fun', resDE.fun], ['GOF', gof], ['Lqq_iv', Lqq_iv],
           ['Lzz_iv', Lzz_iv], ['Crrrr_iv_kPa', Crrrr_iv],
           ['Cqqqq_iv_kPa', Cqqqq_iv], ['Czzzz_iv_kPa', Czzzz_iv],
           ['Wiv_kPa', Wiv], ['CorrMat_det', Corr_Mat_det],
           ['dPdLq', dPdLq], ['PIntercept', PIntercept], 
           ['OD_2_5', OD_2_5], ['OD_5', OD_5], ['OD_7_5', OD_7_5], 
           ['OD_10', OD_10], ['OD_15', OD_15], ['OD_20', OD_20], 
           ['C1', C1], ['C2', C2], ['C3', C3], ['C4', C4], ['C5', C5])


# cd to earlier created folder to store .xlsx file in
os.chdir(z)
# Write created dictionary to .xlsx file
workbook = xlsxwriter.Workbook("%s%s" % (filename, '.xlsx'))
worksheet = workbook.add_worksheet()
row = 0
col = 0
for title, item in results:
    worksheet.write(row, col, title)
    worksheet.write(row+1, col, item)
    col += 1
workbook.close()

''' Save results to MATLAB .mat file '''
# ============================================================================
#mdict = {'resDE': resDE, 'GOF': gof, 'Lqq_iv': Lqq_iv, 'Lzz_iv': Lzz_iv,
#         'Crrrr_iv_kPa': Crrrr_iv, 'Cqqqq_iv_kPa': Cqqqq_iv,
#         'Czzzz_iv_kPa': Czzzz_iv, 'Wiv_kPa': Wiv, 'CorrMat_det': Corr_Mat_det,
#         'dPdLq': dPdLq, 'PIntercept': PIntercept,}
#
#sio.savemat("%s\%s\%s\%s" % (path, filename, 'results'), mdict=mdict,
#            appendmat=True)
# =============================================================================