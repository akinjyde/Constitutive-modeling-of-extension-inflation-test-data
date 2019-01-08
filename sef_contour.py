# -*- coding: utf-8 -*-
''' This module computes the SEF (in vivo), and plots the iso-energy contours'''

# Import needed standard python modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
#import os

def main(x, Lqq_iv, Lzz_iv, model, filename, path):
    
    print('computing the in vivo stored energy density and plotting its contour ...')
    #   kinematics for computing in vivo stored energy density (common to all models)
    Xiv, Yiv = np.meshgrid(Lqq_iv, Lzz_iv)
    Cqqciv = Lqq_iv**2
    Czzciv = Lzz_iv**2
    Crrciv = (1.0/(Lqq_iv*Lzz_iv))**2
    I1civ = Cqqciv + Czzciv + Crrciv
    
    #   kinematics for contour plot of stored energy density (common to all models)
    x = np.linspace(0.9, 1.50, 24)
    y = np.linspace(0.9, 1.50, 24)
    z = 1.0/(x*y)
    #z = np.linspace(0.50, 1.00, 24)
    X, Y = np.meshgrid(x, y)
    Cqqc = X**2
    Czzc = Y**2
    Crrc = z**2
    I1c = Cqqc + Czzc + Crrc
        
    if model == 'NH_2FF':
        
        #   x[0]=ce   x[1]=c1   x[2]=c2   x[3]=alpha        
        I4civ = Cqqciv*(np.sin(x[3]))**2 + Czzciv*(np.cos(x[3]))**2
        Wiv = x[0]*(I1civ - 3.0) + (x[1]/x[2])*(np.exp(x[2]*(I4civ - 1.0)**2) - 1.0)
        print('Wiv [kJ] =', Wiv)
                
        #   Compute stored energy density over the predefined contour area above
        I4c = Cqqc*(np.sin(x[3]))**2 + Czzc*(np.cos(x[3]))**2
        W = x[0]*(I1c - 3.0) + (x[1]/x[2])*(np.exp(x[2]*(I4c - 1.0)**2) - 1.0)
        
    elif model == 'gHY':
        
        #   x[0]=u_T   x[1]=c2   x[2]=E_L   x[3]=u_L   x[4]=c4   x[5]=alpha
        I4civ = Cqqciv*(np.sin(x[5]))**2 + Czzciv*(np.cos(x[5]))**2
        I5civ = (Cqqciv**2)*(np.sin(x[5]))**2 + (Czzciv**2)*(np.cos(x[5]))**2
        A = 0.5*x[0]
        C = 0.5*(x[2] + x[0] - 4.0*x[3])
        B = 0.5*(x[0] - x[3])
        Wiv = ((A/x[1])*(np.exp(x[1]*(I1civ - 3.0)) - 1.0)
               + (C/x[4])*(np.exp(x[4]*(I4civ**0.5 - 1.0)**2) - 1.0)
               + B*(2.0*I4civ - I5civ - 1.0))        
        print('Wiv [kJ] =', Wiv)
                
        #   Compute stored energy density over the predefined contour area above
        I4c = Cqqc*(np.sin(x[5]))**2 + Czzc*(np.cos(x[5]))**2
        I5c = (Cqqc**2)*(np.sin(x[5]))**2 + (Czzc**2)*(np.cos(x[5]))**2
        W = ((A/x[1])*(np.exp(x[1]*(I1c - 3.0)) - 1.0)
               + (C/x[4])*(np.exp(x[4]*(I4c**0.5 - 1.0)**2) - 1.0)
               + B*(2.0*I4c - I5c - 1.0))
        
    elif model == 'gSRM':
        
        #   x[0]=u_T   x[1]=E_L   x[2]=u_L   x[3]=alpha
        I4civ = Cqqciv*(np.sin(x[3]))**2 + Czzciv*(np.cos(x[3]))**2
        I5civ = (Cqqciv**2)*(np.sin(x[3]))**2 + (Czzciv**2)*(np.cos(x[3]))**2
        A = 0.5*x[0]
        C = 0.5*(x[1] + x[0] - 4.0*x[2])
        B = 0.5*(x[0] - x[2])
        Wiv = (A*(I1civ - 3.0)
               + (C/4.0)*(I4civ - 1.0)**2
               + B*(2.0*I4civ - I5civ - 1.0))        
        print('Wiv [kJ] =', Wiv)
                
        #   Compute stored energy density over the predefined contour area above
        I4c = Cqqc*(np.sin(x[3]))**2 + Czzc*(np.cos(x[3]))**2
        I5c = (Cqqc**2)*(np.sin(x[3]))**2 + (Czzc**2)*(np.cos(x[3]))**2
        W = (A*(I1c - 3.0)
               + (C/4.0)*(I4c - 1.0)**2
               + B*(2.0*I4c - I5c - 1.0))   

#    if not os.path.exists(filename):
#      os.makedirs(filename)
#      os.chdir(filename)
      
    #   Make contour plot    
    plt.figure()
    cp = plt.contour(X, Y, W, levels=range(0, 5, 1), colors='k')
    plt.plot(Lqq_iv, Lzz_iv, 'ob', markersize=6, markeredgecolor='b')
    plt.clabel(cp, fmt='%.0f', colors='k', inline=True, fontsize=10)
    plt.gca().yaxis.set_major_formatter(tk.FormatStrFormatter('%.2f'))
    #plt.axes().yaxis.set_major_locator(tk.MultipleLocator(.02))
    plt.gca().xaxis.set_major_formatter(tk.FormatStrFormatter('%.2f'))
    #plt.axes().xaxis.set_major_locator(tk.MultipleLocator(.02))
    plt.title('Contour Plot')
    plt.xlabel(r'$\lambda_\theta$')
    plt.ylabel(r'$\lambda_z$')
    #plt.annotate("Wiv =%sJ" % (Wiv),(1.3,1.45))
    plt.text(1.35, 1.45, r"$W^{}$ = {:2.3f} kJ"
             .format('{iv}', Wiv[0]), color='r')
    plt.text(1.35, 1.41, r"$\lambda_\theta^{}$ = {:4.2f}"
             .format('{iv}', Lqq_iv[0]), color='r')
    plt.text(1.35, 1.37, r"$\lambda_z^{}$ = {:4.2f}"
             .format('{iv}', Lzz_iv[0]), color='r')
    # plt.xlim(1.00, 1.25)
    # plt.ylim(1.00, 1.50)
    plt.savefig("%s\%s\%s\%s" % (path, filename, model, 'contour'), dpi=300, bbox_inches="tight")
    
    return Wiv