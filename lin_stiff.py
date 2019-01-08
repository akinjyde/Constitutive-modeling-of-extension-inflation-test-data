# -*- coding: utf-8 -*-
''' Computes the linearized stifnesses '''

import numpy as np
#import scipy.io as sio

def main(x, Lqq_iv, Lzz_iv, model, filename, path):
    
    print('computing linearized stiffness ...')

    Lrr_iv = 1.0/(Lzz_iv*Lqq_iv)
    np.set_printoptions(precision=4)
    print('Lrr_iv =', Lrr_iv)
    
    Cqq_iv = Lqq_iv**2
    Crr_iv = Lrr_iv**2
    Czz_iv = Lzz_iv**2
    
    I1_iv = Crr_iv + Cqq_iv + Czz_iv
    
    if model == 'NH_2FF':
        
        #      x[0]=ce   x[1]=c1   x[2]=c2   x[3]=alpha
        I4_iv = Cqq_iv*(np.sin(x[3]))**2 + Czz_iv*(np.cos(x[3]))**2
        W1 = 0.5*x[0]
        v = I4_iv - 1.0
        w = x[2]*(v**2)
        
        W4 = 2.0*x[1]*v*np.exp(w)
        
        dW2dLk22 = 4.0*x[1]*x[2]*v*np.exp(w)*(3.0 + 2.0*w)
        
        # Extra stresses in kPa
        Trr_xtra = 2.0*W1*Crr_iv
        Tqq_xtra = 2.0*(W1*Cqq_iv + W4*Cqq_iv*((np.sin(x[3]))**2))
        Tzz_xtra = 2.0*(W1*Czz_iv + W4*Czz_iv*((np.cos(x[3]))**2))
        
        # Linearized stiffnesses in kPa
        Crrrr_iv = 2.0*Trr_xtra        
        Cqqqq_iv = 2.0*Tqq_xtra + 4.0*(Cqq_iv**2)*((np.sin(x[3]))**4)*dW2dLk22
        Czzzz_iv = 2.0*Tzz_xtra + 4.0*(Czz_iv**2)*((np.cos(x[3]))**4)*dW2dLk22
        
        print('Crrrr_iv [kPa] =', Crrrr_iv)
        print('Cqqqq_iv [kPa] =', Cqqqq_iv)
        print('Czzzz_iv [kPa] =', Czzzz_iv)
        
#        mdict = {'Lrr_iv': Lrr_iv, 'Lqq_iv': Lqq_iv, 'Lzz_iv': Lzz_iv,
#                 'Crrrr_iv [kPa]': Crrrr_iv, 'Cqqqq_iv [kPa]': Cqqqq_iv,
#                 'Czzzz_iv [kPa]': Czzzz_iv}
#        
#        sio.savemat("%s\%s\%s\%s" % (path, filename, model, 'LinStiff_results'), mdict=mdict, appendmat=True)   
        
    elif model == 'gHY':
        
        # x[0]=u_T, x[1]=c2, x[2]=E_L, x[3]=u_L, x[4]=c4, and x[5]=alpha
        print('Under development (gHY model)')
        
        # Linearized stiffnesses in kPa
        Crrrr_iv = 0.0        
        Cqqqq_iv = 0.0
        Czzzz_iv = 0.0
        
    elif model == 'gSRM':
        
        # x[0]=u_T, x[1]=E_L, x[2]=u_L, x[3]=alpha
        print('Under development (gSRM model)')
        
        # Linearized stiffnesses in kPa
        Crrrr_iv = 0.0        
        Cqqqq_iv = 0.0
        Czzzz_iv = 0.0
        
    return Crrrr_iv, Cqqqq_iv, Czzzz_iv