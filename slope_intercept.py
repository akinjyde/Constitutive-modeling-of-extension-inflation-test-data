# -*- coding: utf-8 -*-
''' This module computes the slope, and intercept of the P-Lqq curve'''

import numpy as np

def main(x, Lqq_iv, Lzz_iv, ir_iv, or_iv, model, filename):
    
    print('computing Slope and P_Intercept ...')
    
    kPa2Hg = 7.500617

    Lrr_iv = 1.0/(Lzz_iv*Lqq_iv)
        
    Cqq_iv = Lqq_iv**2
    Crr_iv = Lrr_iv**2
    Czz_iv = Lzz_iv**2
    
    I1_iv = Crr_iv + Cqq_iv + Czz_iv
    
    Pe = np.log(or_iv/ir_iv)
    
    if model == 'NH_2FF':
        
        #      x[0]=ce   x[1]=c1   x[2]=c2   x[3]=alpha
        I4_iv = Cqq_iv*(np.sin(x[3]))**2 + Czz_iv*(np.cos(x[3]))**2
        W1 = 0.5*x[0]
        v = I4_iv - 1.0
        w = x[2]*(v**2)
        W4 = 2.0*x[1]*v*np.exp(w)
        dW4dLq = 4.0*x[1]*Lqq_iv*((np.sin(x[3]))**2)*((1.0 + 2.0*w)*np.exp(w))
        
        dPdLq = (2.0*(W1*(2.0*Lqq_iv + ((Lqq_iv**3)*(Lzz_iv**2))**-1)
                 + (2.0*W4*Lqq_iv + Cqq_iv*dW4dLq)*(np.sin(x[3]))**2)*Pe)
        
        PIntercept = 9.0 - (dPdLq*kPa2Hg)*Lqq_iv
        
        print('Slope [kPa] =', dPdLq)
        print('Slope [mmHg] =', dPdLq*kPa2Hg)
        print('P_Intercept [mmHg] =', PIntercept)
        
        
    elif model == 'gHY':
        
        # x[0]=u_T, x[1]=c2, x[2]=E_L, x[3]=u_L, x[4]=c4, and x[5]=alpha
        print('Under development (gHY model)')
        dPdLq = 0
        PIntercept = 0
        
    elif model == 'gSRM':
        
        # x[0]=u_T, x[1]=E_L, x[2]=u_L, x[3]=alpha
        print('Under development (gSRM model)')
        dPdLq = 0
        PIntercept = 0
        
    return dPdLq, PIntercept
    
