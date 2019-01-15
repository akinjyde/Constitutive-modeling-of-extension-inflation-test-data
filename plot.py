# -*- coding: utf-8 -*-
''' Plots figures '''

# Import needed standard python modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tk


def main(x, Lrriv, Lqqiv, Lzziv, iriv, oriv, hiv, ftiv, fiv, OD_iv, P_iv,
         Trriv, Tqqiv, Tzziv, Lrr4a, Lqq4a, Lzz4a, ir4a, or4a, h4a, ft4a,
         f4a, OD_4a, P_4a, Trr4a, Tqq4a, Tzz4a, filename, path, dPdLq,
         PIntercept, Lqq_iv):
    
    print('making plots ...')
    # Unit conversion parameters
    kPa2Hg = 7.500617
    kN2mN = 1.0e6
    N2mN = 1.0e3
    Pa2kPa = 1.0e-3
    Pa2Hg = 0.00750062
    
   
    Trriv = Trriv*Pa2kPa
    Tqqiv = Tqqiv*Pa2kPa
    Tzziv = Tzziv*Pa2kPa
    
    Trr4a = Trr4a*Pa2kPa
    Tqq4a = Tqq4a*Pa2kPa
    Tzz4a = Tzz4a*Pa2kPa
    # =========================================================================
    
    Crriv = Lrriv**2
    Cqqiv = Lqqiv**2
    Czziv = Lzziv**2
    
    Crr4a = Lrr4a**2
    Cqq4a = Lqq4a**2
    Czz4a = Lzz4a**2
    
  
    # Import model-dependent module, and functions
    from theory_nh2ff import Pressure_th, TForce_th
    
    # Compute theoretical for in vivo with optimized parameters
    
    P_th_iv = Pressure_th(x, Crriv, Cqqiv, Czziv, iriv, oriv)   # kPa        
    ft_th_iv = TForce_th(x, Crriv, Cqqiv, Czziv, iriv, oriv)    # kN                
    fp_th_iv = P_th_iv*np.pi*iriv**2                            # kN
    f_th_iv = fp_th_iv + ft_th_iv                               # kN
    
    Trrivth = (-P_th_iv*iriv)/(oriv + iriv)                     # kPa
    Tqqivth = (iriv*P_th_iv)/hiv                                # kPa
    Tzzivth = f_th_iv/(np.pi*(oriv**2 - iriv**2))               # kPa
    # Addition of shear stress made on February 3rd, 2018
    I4 = Cqqiv*((np.sin(x[3]))**2) + Czziv*((np.cos(x[3]))**2)
    v = I4 - 1.0
    w = x[2]*(v**2)                                             # Unitless        
    W4 = 2.0*x[1]*v*np.exp(w)                                   # kPa
    Tqzivth = 2.0*W4*Lqqiv*Lzziv*np.sin(x[3])*np.cos(x[3])      # kPa ?

    # Compute theoretical for 4% ABOVE with optimized parameters
   
    P_th_4a = Pressure_th(x, Crr4a, Cqq4a, Czz4a, ir4a, or4a)   # kPa       
    ft_th_4a = TForce_th(x, Crr4a, Cqq4a, Czz4a, ir4a, or4a)    # kN
    fp_th_4a = P_th_4a*np.pi*ir4a**2                            # kN
    f_th_4a = fp_th_4a + ft_th_4a                               # kN
    
    Trr4ath = (-P_th_4a*ir4a)/(or4a + ir4a)                     # kPa
    Tqq4ath = (ir4a*P_th_4a)/h4a                                # kPa
    Tzz4ath = f_th_4a/(np.pi*(or4a**2 - ir4a**2))               # kPa
        
    
    # =========================================================================
    ''' Pressure vs outer diameter plot (in vivo and 4% above) '''
    plt.figure()
#    # +4%
#    plt.plot(OD_4a, P_4a*Pa2Hg, 's', color='k', markersize=6,
#             fillstyle='none', label='exp (iv+4%)')
#    plt.plot(OD_4a, P_th_4a*kPa2Hg, color='k')
    # in vivo
    plt.plot(OD_iv, P_iv*Pa2Hg, 'o', markersize=6,
             fillstyle='none', color='k', label='exp (iv)')
    plt.plot(OD_iv, P_th_iv*kPa2Hg, color='k', label='model')
    # Turn on grid for slope estimation
    plt.grid(b=True, which='both', axis='both', color='b', linestyle='--',
             linewidth=0.25)
    
    plt.xlabel('Outer Diameter ($\mu$m)', fontsize=12)
    plt.ylabel('Pressure (mmHg)', fontsize=12)
    plt.title(filename, fontsize=12)
    plt.legend(loc='upper center')
    plt.savefig("%s\%s\%s" % (path, filename, 'P_OD'), dpi=300,
                bbox_inches="tight")
    # ========================================================================
    plt.figure()
    plt.plot(Lqqiv, P_iv*Pa2Hg, 'o', markersize=6,
             fillstyle='none', color='k', label='exp (iv)')
    plt.plot(Lqqiv, P_th_iv*kPa2Hg, color='k', label='model')
    Lqq_tangent = np.array([Lqq_iv-0.005, Lqq_iv+0.005])
    plt.plot(Lqq_tangent, (dPdLq*kPa2Hg)*Lqq_tangent + PIntercept, color='r',
             lw=1.5)
    plt.axhline(y=9.0, linewidth=0.5, linestyle='--', color='r')
    plt.axvline(x=Lqq_iv, linewidth=0.5, linestyle='--', color='r')
    
    plt.xlabel('Outer Diameter ($\mu$m)', fontsize=12)
    plt.ylabel('Pressure (mmHg)', fontsize=12)
    plt.title(filename, fontsize=12)
    plt.legend(loc='upper center')
    plt.savefig("%s\%s\%s" % (path, filename, 'P_Lqq_tangent'),
                dpi=300, bbox_inches="tight")
    # =========================================================================
    ''' Total Axial force vs pressure (in vivo and 4% above) '''
    plt.figure()
#    # +4%
#    plt.plot(P_4a*Pa2Hg, f4a*N2mN, 's', color='k', markersize=6,
#             fillstyle='none', label='exp (iv+4%)')
#    plt.plot(P_4a*Pa2Hg, f_th_4a*kN2mN, color='k')
    # in vivo
    plt.plot(P_iv*Pa2Hg, fiv*N2mN, 'o', markersize=6,
             fillstyle='none', color='k', label='exp (iv)')
    plt.plot(P_iv*Pa2Hg, f_th_iv*kN2mN, color='k', label='model')
    
    plt.xlabel('Pressure (mmHg)', fontsize=12)
    plt.ylabel('Total axial Force (mN)', fontsize=12)
    plt.title(filename, fontsize=12)
    plt.legend(loc='upper center')
    plt.savefig("%s\%s\%s" % (path, filename, 'F_P'), dpi=300,
                bbox_inches="tight")
    # =========================================================================
    ''' Transducer-measured force vs pressure (in vivo and 4% above) '''
    plt.figure()
#    # +4%
#    plt.plot(P_4a*Pa2Hg, ft4a*N2mN, 's', color='k', markersize=6,
#             fillstyle='none', label='exp (iv+4%)')
#    plt.plot(P_4a*Pa2Hg, ft_th_4a*kN2mN, color='k')
    # in vivo
    plt.plot(P_iv*Pa2Hg, ftiv*N2mN, 'o', markersize=6,
             fillstyle='none', color='k', label='exp (iv)')
    plt.plot(P_iv*Pa2Hg, ft_th_iv*kN2mN, color='k', label='model')
    
    plt.xlabel('Pressure (mmHg)', fontsize=12)
    plt.ylabel('Axial Force (mN)', fontsize=12)
    plt.title(filename, fontsize=12)
    plt.legend(loc='upper center')
    plt.savefig("%s\%s\%s" % (path, filename, 'FT_P'), dpi=300,
                bbox_inches="tight")
    # =========================================================================
    ''' Tqq vs Lqq plot (in vivo and 4% above) '''
    plt.figure()
#    # +4%
#    plt.plot(Lqq4a, Tqq4a, 's', color='k', markersize=6, fillstyle='none',
#             label='exp (iv+4%)')
#    plt.plot(Lqq4a, Tqq4ath, color='k')
    # in vivo
    plt.plot(Lqqiv, Tqqiv, 'o', color='k', markersize=6, fillstyle='none',
             label='exp (iv)')
    plt.plot(Lqqiv, Tqqivth, color='k', label='model')
    
    plt.gca().xaxis.set_major_formatter(tk.FormatStrFormatter('%.2f'))
    plt.axes().xaxis.set_major_locator(tk.MultipleLocator(.05))
    #plt.xlim(0.85, 1.16)
    plt.xlabel(r'$\lambda_\theta$', fontsize=12)
    plt.ylabel(r'$\sigma_{\theta\theta}$ (kPa)', fontsize=12)
    plt.title(filename, fontsize=12)
    plt.legend(loc='upper center')
    plt.savefig("%s\%s\%s" % (path, filename, 'Tqq_Lqq'), dpi=300,
                bbox_inches="tight")
    # =========================================================================
    ''' Tzz vs Lqq plot (in vivo and 4% above) '''
    plt.figure()
#    # +4%
#    plt.plot(Lqq4a, Tzz4a, 's', color='k', markersize=6, fillstyle='none',
#             label='exp (iv+4%)')
#    plt.plot(Lqq4a, Tzz4ath, color='k')
    # in vivo
    plt.plot(Lqqiv, Tzziv, 'o', color='k', markersize=6, fillstyle='none',
             label='exp (iv)')
    plt.plot(Lqqiv, Tzzivth, color='k', label='model')
    
    plt.gca().xaxis.set_major_formatter(tk.FormatStrFormatter('%.2f'))
    plt.axes().xaxis.set_major_locator(tk.MultipleLocator(.02))
    plt.xlabel(r'$\lambda_\theta$', fontsize=12)
    plt.ylabel(r'$\sigma_{zz}$ (kPa)', fontsize=12)
    plt.title(filename, fontsize=12)
    plt.legend(loc='upper center')
    plt.savefig("%s\%s\%s" % (path, filename, 'Tzz_Lqq'), dpi=300,
                bbox_inches="tight")
    # =========================================================================
    ''' Tzz vs Lzz plot (in vivo and 4% above) '''
    #plt.figure()
    #plt.plot(Lzz4a, Tzz4a, 's', color='k', markersize=6, fillstyle='none',
    #         label='exp (iv+4%)')
    #plt.plot(Lzz4a, Tzz4ath, color='k')
    #
    #plt.plot(Lzziv, Tzziv, 'o', color='k', markersize=6, fillstyle='none',
    #         label='exp (iv)')
    #plt.plot(Lzziv, Tzzivth, color='k', label='model')
    #
    #plt.gca().xaxis.set_major_formatter(tk.FormatStrFormatter('%.2f'))
    #plt.axes().xaxis.set_major_locator(tk.MultipleLocator(.02))
    #plt.xlabel(r'$\lambda_z$', fontsize=12)
    #plt.ylabel(r'$\sigma_{zz}$ (kPa)', fontsize=12)
    #plt.title(filename, fontsize=12)
    #plt.legend(loc='upper center')
    #plt.savefig("%s\%s\%s\%s" % (path, filename, model, 'Tzz_Lzz'), dpi=300, bbox_inches="tight")
    # =========================================================================
    ''' Trr vs Lrr plot (in vivo and 4% above) '''
    plt.figure()
#    # +4%
#    plt.plot(Lrr4a, Trr4a, 's', color='k', markersize=6, fillstyle='none',
#             label='exp (iv+4%)')
#    plt.plot(Lrr4a, Trr4ath, color='k')
    # in vivo
    plt.plot(Lrriv, Trriv, 'o', color='k', markersize=6, fillstyle='none',
             label='exp (iv)')
    plt.plot(Lrriv, Trrivth, color='k', label='model')
    
    plt.gca().xaxis.set_major_formatter(tk.FormatStrFormatter('%.2f'))
    plt.axes().xaxis.set_major_locator(tk.MultipleLocator(.02))
    plt.xlabel(r'$\lambda_r$', fontsize=12)
    plt.ylabel(r'$\sigma_{rr}$ (kPa)', fontsize=12)
    plt.title(filename, fontsize=12)
    plt.legend(loc='lower center')
    plt.savefig("%s\%s\%s" % (path, filename, 'Trr_Lrr'), dpi=300,
                bbox_inches="tight")
    
    return (Trrivth, Tqqivth, Tzzivth, Tqzivth, P_th_iv, ft_th_iv, f_th_iv, Trr4ath,
            Tqq4ath, Tzz4ath, P_th_4a, ft_th_4a, f_th_4a)