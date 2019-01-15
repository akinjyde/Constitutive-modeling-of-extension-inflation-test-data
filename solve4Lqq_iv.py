# -*- coding: utf-8 -*-
''' This module solves for the in vivo circumferential stretch at 17mmHg '''

# Import needed standard python modules
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import matplotlib.ticker as tk


def main(x, P_iv, Lqqiv, Lzziv, iriv, oriv, filename, path):
    
    print('solving for the in vivo circumferential stretch ...')
    
    # 1 kPa = 7.500617 mmHg, convert P_iv from kPa to mmHg
    kPa2Hg = 7.500617
    Pa2Hg = 0.00750062
    # =============================================================================
    # Interpolate the in vivo circumferential stretch at 17 mmHg
    Lqq_iv_est = np.interp(9.0, np.sort(np.ravel(P_iv*Pa2Hg), axis=None),
                           np.sort(np.ravel(Lqqiv), axis=None))
    # Interpolate the in vivo inner radius at 17 mmHg
    ir_iv = np.interp(9.0, np.sort(np.ravel(P_iv*Pa2Hg), axis=None),
                      np.sort(np.ravel(iriv), axis=None))
    # Interpolate the in vivo outer radius at 17 mmHg
    or_iv = np.interp(9.0, np.sort(np.ravel(P_iv*Pa2Hg), axis=None),
                      np.sort(np.ravel(oriv), axis=None))
    # =============================================================================
    # Find the in vivo circumferential stretch
    print('Lqq_iv_est =', Lqq_iv_est)
    Lzz_iv = Lzziv[1]
    np.set_printoptions(precision=4)
    print('Lzz_iv =', Lzz_iv)
    
    z0 = np.array([Lqq_iv_est])
      
        
    def P_Lqiv(z, x, Lzz_iv, ir_iv, or_iv, Pa2Hg):
        
        Lqq_iv = z[0]
        Cqq_iv = Lqq_iv**2
        Crr_iv = (1.0/(Lqq_iv*Lzz_iv))**2
        Czz_iv = Lzz_iv**2
        
        
        #      x[0]=ce   x[1]=c1   x[2]=c2   x[3]=alpha
        I4_iv = Cqq_iv*(np.sin(x[3]))**2 + Czz_iv*(np.cos(x[3]))**2
        W1 = 0.5*(x[0]*Pa2Hg)
        W4 = 2.0*(x[1]*Pa2Hg)*(I4_iv - 1.0)*np.exp(x[2]*(I4_iv - 1.0)**2)
        P_th = 2.0*(W1*(Cqq_iv - Crr_iv) + W4*Cqq_iv*(np.sin(x[3]))**2)*np.log(or_iv/ir_iv)
            
        fun = P_th - 17.0
        return fun
        
    Lqq_iv = fsolve(P_Lqiv, z0, args=(x, Lzz_iv, ir_iv, or_iv, Pa2Hg))
    print('Lqq_iv =', Lqq_iv)
    
    # =============================================================================
    # Piv vs Lqq plot to display solved Lqq
    plt.figure()
    plt.plot(Lqqiv, P_iv*Pa2Hg, 'o:', markersize=6, fillstyle='none', color='b',
             label=r'$P^{iv} (exp)$')
    plt.text((Lqq_iv - 0.04), 9.5, '9 mmHg', color='r')
    plt.text((Lqq_iv - 0.02), 0, r"$\lambda_\theta^{}$ = {:5.3f}"
             .format('{iv}', Lqq_iv[0, ]), color='r')
    plt.axhline(y=9, linewidth=0.5, linestyle='--', color='r')
    plt.axvline(x=Lqq_iv, linewidth=0.5, linestyle='--', color='r')
    
    plt.gca().xaxis.set_major_formatter(tk.FormatStrFormatter('%.2f'))
    plt.axes().xaxis.set_major_locator(tk.MultipleLocator(.02))
    plt.xlabel(r'$\lambda_\theta$', fontsize=12)
    plt.ylabel("P_iv (mmHg)", fontsize=12)
    plt.title(filename, fontsize=12)
    plt.legend(loc='upper center')
    plt.savefig("%s\%s\%s" % (path, filename, 'Piv_Lqiv'), dpi=300, bbox_inches="tight")
    # =============================================================================
    return Lqq_iv, Lzz_iv, ir_iv, or_iv