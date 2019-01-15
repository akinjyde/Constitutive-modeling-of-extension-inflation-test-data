# -*- coding: utf-8 -*-
''' This module compute the correlation matrix and the local sensitivity'''

# Import needed standard python modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
import seaborn as sns

def main(x, Crr, Cqq, Czz, Or, Ir, h, Tqqivth, Tzzivth, filename, path):
    
    print('Computing correlation matrix and sensitivity analysis')
        
    Pe = np.log(Or/Ir)
    Le = 0.5*np.pi*(Or**2 - Ir**2)
    me = 1.0/(np.pi*(Or**2 - Ir**2))
    
    # Correlation matrix
    
        
    #   x[0]=ce   x[1]=c1   x[2]=c2   x[3]=alpha
    I4 = Cqq*((np.sin(x[3]))**2) + Czz*((np.cos(x[3]))**2)
    v = I4 - 1.0
    w = x[2]*(v**2)
    
    W4 = 2.0*x[1]*v*np.exp(w)
    
    dPdW4 = 2.0*Cqq*((np.sin(x[3]))**2)*Pe                  #   dPdW4        
    dW4dc1 = 2.0*v*np.exp(w)                                #   dW4dx[1]
    dW4dc2 = 2.0*x[1]*(v**3)*np.exp(w)                      #   dW4dx[2]
    dW4da = 2.0*x[1]*(Cqq - Czz)*np.sin(2*x[3])*np.exp(w)*(1.0 + 2.0*w)  # dW4dx[3]
    
    #   change in intraluminal pressure wrt parameters
    dPdc = (Cqq - Crr)*Pe                                   #   dPdx[0]
    dPdc1 = dPdW4*dW4dc1                                    #   dPdx[1]
    dPdc2 = dPdW4*dW4dc2                                    #   dPdx[2]
    dPda = 2.0*Cqq*(W4*np.sin(2*x[3]) + ((np.sin(x[3]))**2)*dW4da)*Pe  # dPdx[3]
            
    be = 2.0*Czz*((np.cos(x[3]))**2) - Cqq*((np.sin(x[3]))**2)   
    
    #   change in transducer-measured force wrt parameters
    dLdc = (2.0*Czz - Crr - Cqq)*Le
    dLdc1 = 2.0*dW4dc1*be*Le
    dLdc2 = 2.0*dW4dc2*be*Le
    dLda = 2.0*(dW4da*be - W4*np.sin(2*x[3])*(2.0*Czz + Cqq))*Le
            
    #   change in circumferential stress wrt parameters
    dTqdc = dPdc*Ir/h
    dTqdc1 = dPdc1*Ir/h
    dTqdc2 = dPdc2*Ir/h
    dTqda = dPda*Ir/h
    
    #   change in axial stress wrt model parameters
    dTzdc = (np.pi*(Ir**2)*dPdc + dLdc)*me
    dTzdc1 = (np.pi*(Ir**2)*dPdc1 + dLdc1)*me
    dTzdc2 = (np.pi*(Ir**2)*dPdc2 + dLdc2)*me
    dTzda = (np.pi*(Ir**2)*dPda + dLda)*me
    
    J_0 = np.concatenate((dTqdc, dTzdc), axis=0)
    J_1 = np.concatenate((dTqdc1, dTzdc1), axis=0)
    J_2 = np.concatenate((dTqdc2, dTzdc2), axis=0)
    J_3 = np.concatenate((dTqda, dTzda), axis=0)
    
    #print('J_0.shape =', J_0.shape)
    J = np.empty(((Czz.size)*2, x.size))
    #print('J.shape =', J.shape)
    J[:, 0] = np.ravel(J_0)
    J[:, 1] = np.ravel(J_1)
    J[:, 2] = np.ravel(J_2)
    J[:, 3] = np.ravel(J_3)
    J = np.nan_to_num(J)
    
    Corr_Mat = np.corrcoef(J, rowvar=False)
    Corr_Mat_det = np.linalg.det(Corr_Mat)
    print('Corr Mat =', Corr_Mat)
    np.set_printoptions(precision=4)
    print('Det Corr Mat =', Corr_Mat_det)
    if Corr_Mat_det < 1e-4:
        print('Overparameterized!')
    
    plt.figure()
    ax = sns.heatmap(Corr_Mat, vmin=-1.0, vmax=1.0, cmap=plt.cm.RdBu, center=0.0,
                     robust=False, annot=True, fmt='.3g', linewidths=.5,
                     xticklabels=[r'$\mu$', r'$c_1$', r'$c_2$', r'$\alpha$'],
                     yticklabels=[r'$\mu$', r'$c_1$', r'$c_2$', r'$\alpha$'])
    plt.title(filename, fontsize=12)
    plt.savefig("%s\%s\%s" % (path, filename, 'Corr'), dpi=300, bbox_inches="tight") 
    # =============================================================================
    # Plot Sensitivity Analaysis with subplot
    plt.figure()
    plt.subplot(2,1,1)
    #   Circumferential Stress
    Lqq = np.sqrt(Cqq)
    plt.plot(Lqq, dTqdc*(x[0]/Tqqivth), color='k', label=r'$\mu$')
    plt.plot(Lqq, dTqdc1*(x[1]/Tqqivth), color='r', label=r'$c_1$')
    plt.plot(Lqq, dTqdc2*(x[2]/Tqqivth), color='b', label=r'$c_2$')
    plt.plot(Lqq, dTqda*(x[3]/Tqqivth), color='g', label=r'$\alpha$')
    plt.gca().yaxis.set_major_formatter(tk.FormatStrFormatter('%.1f'))
    #plt.axes().yaxis.set_major_locator(tk.MultipleLocator(.02))
    #plt.xlabel(r'$\lambda_\theta$', fontsize=12)
    #plt.legend(loc=9, bbox_to_anchor=(0.4, -0.5), ncol=4)
    #plt.legend(bbox_to_anchor=(1.04, -0.01), loc=2, borderaxespad=0.)
    plt.ylabel(r'$S^{\sigma_\theta}$', fontsize=12)
    plt.title(filename, fontsize=8)
    plt.legend(loc='best', ncol=4)
    
    plt.subplot(2,1,2)
    #   Axial Stress
    plt.plot(Lqq, dTzdc*(x[0]/Tzzivth), color='k', label=r'$\mu$')
    plt.plot(Lqq, dTzdc1*(x[1]/Tzzivth), color='r', label=r'$c_1$')
    plt.plot(Lqq, dTzdc2*(x[2]/Tzzivth), color='b', label=r'$c_2$')
    plt.plot(Lqq, dTzda*(x[3]/Tzzivth), color='g', label=r'$\alpha$')
    plt.gca().yaxis.set_major_formatter(tk.FormatStrFormatter('%.1f'))
    #plt.axes().xaxis.set_major_locator(tk.MultipleLocator(.02))
    plt.xlabel(r'$\lambda_\theta$', fontsize=12)
    plt.ylabel(r'$S^{\sigma_z}$', fontsize=12)
    #plt.legend(loc='best')
    plt.savefig("%s\%s\%s" % (path, filename, 'Sens'), dpi=300, bbox_inches="tight")
        
    return Corr_Mat_det
    
    