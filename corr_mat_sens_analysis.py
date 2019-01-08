# -*- coding: utf-8 -*-
''' This module compute the correlation matrix and the local sensitivity'''

# Import needed standard python modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
import seaborn as sns

def main(x, Crr, Cqq, Czz, Or, Ir, h, Tqqivth, Tzzivth, model, filename, path):
    
    print('Computing correlation matrix and sensitivity analysis')
        
    Pe = np.log(Or/Ir)
    Le = 0.5*np.pi*(Or**2 - Ir**2)
    me = 1.0/(np.pi*(Or**2 - Ir**2))
    
    # Correlation matrix
    if model == 'NH_2FF':
        
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
        plt.savefig("%s\%s\%s\%s" % (path, filename, model, 'Corr'), dpi=300, bbox_inches="tight") 
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
        plt.savefig("%s\%s\%s\%s" % (path, filename, model, 'Sens'), dpi=300, bbox_inches="tight")
        
    elif model == 'gHY':
        
        #   x[0]=u_T   x[1]=c2   x[2]=E_L   x[3]=u_L   x[4]=c4   x[5]=alpha
        
        I1 = Crr + Cqq + Czz
        I4 = Cqq*(np.sin(x[5]))**2 + Czz*(np.cos(x[5]))**2
        u = I4**-0.5
        duda = -0.5*(u**3)*np.sin(2*x[5])*(Cqq - Czz)
        v = I4**0.5 - 1.0
        dvda = 0.5*u*np.sin(2*x[5])*(Cqq - Czz)
        w2 = x[1]*(I1 - 3.0)
        w4 = x[4]*(v**2)
        
        A = 0.5*x[0]
        C = 0.5*(x[2] + x[0] - 4.0*x[3])
        B = 0.5*(x[0] - x[3])
        W1 = 0.5*A*np.exp(w2)
        W4 = C*(I4**-0.5)*v*np.exp(w4) + 2.0*B
        W5 = -B
        
        dW1duT = 0.25*x[1]*np.exp(w2)                                #   dW1dx[1]
        dW1dc2 = 0.5*A*np.exp(w2)*(1.0 + w2)                      #   dW4dx[2]
        dW4duT = 0.5*v*u*np.exp(w4) + 1.0
        dW4dEL = 0.5*v*u*np.exp(w4)
        dW4duL = -(2*v*u*np.exp(w4) + 1.0)
        dW4dc4 = C*(v**3)*u*np.exp(w4)
        dW4da = C*(u*dvda*np.exp(w4)*(2.0*w4 + 1.0) + v*np.exp(w4)*duda)
        dW5duT = -0.5
        dW5duL = 0.5
        
        #   change in intraluminal pressure wrt parameters
        dPduT = 2.0*(dW1duT*(Cqq - Crr) + Cqq*(np.sin(x[5]))**2*(dW4duT + 2.0*dW5duT*Cqq))*Pe                  #   dPdW4        
        dPdc2 = 2.0*(dW1dc2*(Cqq - Crr))*Pe
        dPdEL = 2.0*(Cqq*((np.sin(x[5]))**2)*dW4dEL)*Pe                                    #   dPdx[2]
        dPduL = 2.0*(Cqq*((np.sin(x[5]))**2)*(dW4duL + 2.0*dW5duL*Cqq))*Pe
        dPdc4 = 2.0*(Cqq*((np.sin(x[5]))**2)*dW4dc4)*Pe
        dPda = 2.0*(Cqq*((np.sin(x[5]))**2)*dW4da)*Pe
                
        be = 2.0*Czz*((np.cos(x[5]))**2) - Cqq*((np.sin(x[5]))**2)   
        bf = 2.0*(Czz**2)*((np.cos(x[5]))**2) - (Cqq**2)*((np.sin(x[5]))**2)
                
        #   change in transducer-measured force wrt parameters
        dLduT = (dW1duT*(2.0*Czz - Crr - Cqq) + dW4duT*be + 2.0*dW5duT*bf)*Le
        dLdc2 = 2.0*(dW1dc2*be)*Le
        dLdEL = 2.0*(dW4dEL*be)*Le
        dLduL = 2.0*(dW4duL*be + 2.0*dW5duL*bf)*Le
        dLdc4 = 2.0*(dW4dc4*be)*Le
        dLda = 2.0*(dW4da*be - W4*np.sin(2.0*x[5])*(2.0*Czz + Cqq))*Le
                
        #   change in circumferential stress wrt parameters
        dTqduT = dPduT*Ir/h
        dTqdc2 = dPdc2*Ir/h
        dTqdEL = dPdEL*Ir/h
        dTqduL = dPduL*Ir/h
        dTqdc4 = dPdc4*Ir/h
        dTqda = dPda*Ir/h
        
        #   change in axial stress wrt model parameters
        dTzduT = (np.pi*(Ir**2)*dPduT + dLduT)*me
        dTzdc2 = (np.pi*(Ir**2)*dPdc2 + dLdc2)*me
        dTzdEL = (np.pi*(Ir**2)*dPdEL + dLdEL)*me
        dTzduL = (np.pi*(Ir**2)*dPduL + dLduL)*me
        dTzdc4 = (np.pi*(Ir**2)*dPdc4 + dLdc4)*me
        dTzda = (np.pi*(Ir**2)*dPda + dLda)*me
        
        J_0 = np.concatenate((dTqduT, dTzduT), axis=0)
        J_1 = np.concatenate((dTqdc2, dTzdc2), axis=0)
        J_2 = np.concatenate((dTqdEL, dTzdEL), axis=0)
        J_3 = np.concatenate((dTqduL, dTzduL), axis=0)
        J_4 = np.concatenate((dTqdc4, dTzdc4), axis=0)
        J_5 = np.concatenate((dTqda, dTzda), axis=0)
        
        #print('J_0.shape =', J_0.shape)
        J = np.empty(((Czz.size)*2, x.size))
        #print('J.shape =', J.shape)
        J[:, 0] = np.ravel(J_0)
        J[:, 1] = np.ravel(J_1)
        J[:, 2] = np.ravel(J_2)
        J[:, 3] = np.ravel(J_3)
        J[:, 4] = np.ravel(J_4)
        J[:, 5] = np.ravel(J_5)
        
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
                         xticklabels=[r'$\mu_T$', r'$c_2$', r'$E_L$', r'$\mu_L$', r'$c_4$', r'$\alpha$'],
                         yticklabels=[r'$\mu_T$', r'$c_2$', r'$E_L$', r'$\mu_L$', r'$c_4$', r'$\alpha$'])
        plt.title(filename, fontsize=12)
        plt.savefig("%s\%s\%s\%s" % (path, filename, model, 'Corr'), dpi=300, bbox_inches="tight") 
        # =============================================================================
        # Plot Sensitivity Analaysis with subplot
        plt.figure()
        plt.subplot(2,1,1)
        Lqq = np.sqrt(Cqq)
        plt.plot(Lqq, dTqduT*(x[0]/Tqqivth), color='k', label=r'$\mu_T$')
        plt.plot(Lqq, dTqdc2*(x[1]/Tqqivth), color='r', label=r'$c_2$')
        plt.plot(Lqq, dTqdEL*(x[2]/Tqqivth), color='b', label=r'$E_L$')
        plt.plot(Lqq, dTqduL*(x[3]/Tqqivth), color='c', label=r'$\mu_L$')
        plt.plot(Lqq, dTqdc4*(x[4]/Tqqivth), color='m', label=r'$c_4$')
        plt.plot(Lqq, dTqda*(x[5]/Tqqivth), color='g', label=r'$\alpha$')
        plt.gca().yaxis.set_major_formatter(tk.FormatStrFormatter('%.1f'))
        #plt.axes().yaxis.set_major_locator(tk.MultipleLocator(.02))
        #plt.xlabel(r'$\lambda_\theta$', fontsize=12)
        #plt.legend(loc=9, bbox_to_anchor=(0.4, -0.5), ncol=4)
        #plt.legend(bbox_to_anchor=(1.04, -0.01), loc=2, borderaxespad=0.)
        plt.ylabel(r'$S^{\sigma_\theta}$', fontsize=12)
        plt.title(filename, fontsize=8)
        plt.legend(loc='best', ncol=4)
        
        plt.subplot(2,1,2)
        plt.plot(Lqq, dTzduT*(x[0]/Tzzivth), color='k')
        plt.plot(Lqq, dTzdc2*(x[1]/Tzzivth), color='r')
        plt.plot(Lqq, dTzdEL*(x[2]/Tzzivth), color='b')
        plt.plot(Lqq, dTzduL*(x[3]/Tzzivth), color='c')
        plt.plot(Lqq, dTzdc4*(x[4]/Tzzivth), color='m')
        plt.plot(Lqq, dTzda*(x[5]/Tzzivth), color='g')
        plt.gca().yaxis.set_major_formatter(tk.FormatStrFormatter('%.1f'))
        #plt.axes().xaxis.set_major_locator(tk.MultipleLocator(.02))
        plt.xlabel(r'$\lambda_\theta$', fontsize=12)
        plt.ylabel(r'$S^{\sigma_z}$', fontsize=12)
        #plt.legend(loc='best')
        plt.savefig("%s\%s\%s\%s" % (path, filename, model, 'Sens'), dpi=300, bbox_inches="tight")
        
    elif model == 'gSRM':
        
        #   x[0]=u_T   x[1]=E_L   x[2]=u_L   x[3]=alpha
        
        I1 = Crr + Cqq + Czz
        I4 = Cqq*(np.sin(x[3]))**2 + Czz*(np.cos(x[3]))**2
        v = I4 - 1.0
               
        A = 0.5*x[0]
        C = 0.5*(x[1] + x[0] - 4.0*x[2])
        B = 0.5*(x[0] - x[2])
        W1 = A
        W4 = 0.5*C*(I4 - 1.0) + 2.0*B
        W5 = -B
        
        dW1duT = 0.5                                #   dW1dx[1]
        dW4duT = 0.25*v + 1.0
        dW4dEL = 0.25*v
        dW4duL = I4
        dW4da = 0.5*C*np.sin(2*x[3])*(Cqq - Czz)
        dW5duT = -0.5
        dW5duL = 0.5
        
        #   change in intraluminal pressure wrt parameters
        dPduT = 2.0*(dW1duT*(Cqq - Crr) + Cqq*(np.sin(x[3]))**2*(dW4duT + 2.0*dW5duT*Cqq))*Pe                  #   dPdW4        
        dPdEL = 2.0*(Cqq*((np.sin(x[3]))**2)*dW4dEL)*Pe                                    #   dPdx[2]
        dPduL = 2.0*(Cqq*((np.sin(x[3]))**2)*(dW4duL + 2.0*dW5duL*Cqq))*Pe
        dPda = 2.0*(Cqq*((np.sin(x[3]))**2)*dW4da)*Pe
                
        be = 2.0*Czz*((np.cos(x[3]))**2) - Cqq*((np.sin(x[3]))**2)   
        bf = 2.0*(Czz**2)*((np.cos(x[3]))**2) - (Cqq**2)*((np.sin(x[3]))**2)
                
        #   change in transducer-measured force wrt parameters
        dLduT = (dW1duT*(2.0*Czz - Crr - Cqq) + dW4duT*be + 2.0*dW5duT*bf)*Le
        dLdEL = 2.0*(dW4dEL*be)*Le
        dLduL = 2.0*(dW4duL*be + 2.0*dW5duL*bf)*Le
        dLda = 2.0*(dW4da*be - W4*np.sin(2.0*x[3])*(2.0*Czz + Cqq))*Le
                
        #   change in circumferential stress wrt parameters
        dTqduT = dPduT*Ir/h
        dTqdEL = dPdEL*Ir/h
        dTqduL = dPduL*Ir/h
        dTqda = dPda*Ir/h
        
        #   change in axial stress wrt model parameters
        dTzduT = (np.pi*(Ir**2)*dPduT + dLduT)*me
        dTzdEL = (np.pi*(Ir**2)*dPdEL + dLdEL)*me
        dTzduL = (np.pi*(Ir**2)*dPduL + dLduL)*me
        dTzda = (np.pi*(Ir**2)*dPda + dLda)*me
        
        J_0 = np.concatenate((dTqduT, dTzduT), axis=0)
        J_1 = np.concatenate((dTqdEL, dTzdEL), axis=0)
        J_2 = np.concatenate((dTqduL, dTzduL), axis=0)
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
                         xticklabels=[r'$\mu_T$', r'$E_L$', r'$\mu_L$', r'$\alpha$'],
                         yticklabels=[r'$\mu_T$', r'$E_L$', r'$\mu_L$', r'$\alpha$'])
        plt.title(filename, fontsize=12)
        plt.savefig("%s\%s\%s\%s" % (path, filename, model, 'Corr'), dpi=300, bbox_inches="tight") 
        # =============================================================================
        # Plot Sensitivity Analaysis with subplot
        plt.figure()
        plt.subplot(2,1,1)
        Lqq = np.sqrt(Cqq)
        plt.plot(Lqq, dTqduT*(x[0]/Tqqivth), color='k', label=r'$\mu_T$')
        plt.plot(Lqq, dTqdEL*(x[1]/Tqqivth), color='b', label=r'$E_L$')
        plt.plot(Lqq, dTqduL*(x[2]/Tqqivth), color='c', label=r'$\mu_L$')
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
        plt.plot(Lqq, dTzduT*(x[0]/Tzzivth), color='k')
        plt.plot(Lqq, dTzdEL*(x[1]/Tzzivth), color='b')
        plt.plot(Lqq, dTzduL*(x[2]/Tzzivth), color='c')
        plt.plot(Lqq, dTzda*(x[3]/Tzzivth), color='g')
        plt.gca().yaxis.set_major_formatter(tk.FormatStrFormatter('%.1f'))
        #plt.axes().xaxis.set_major_locator(tk.MultipleLocator(.02))
        plt.xlabel(r'$\lambda_\theta$', fontsize=12)
        plt.ylabel(r'$S^{\sigma_z}$', fontsize=12)
        #plt.legend(loc='best')
        plt.savefig("%s\%s\%s\%s" % (path, filename, model, 'Sens'), dpi=300, bbox_inches="tight")
        
    return Corr_Mat_det
    
    