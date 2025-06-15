import numpy as np
import csv
import pandas as pd
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import h5py
import math
import array
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import ticker
from matplotlib import rcParams
import matplotlib as mpl
mpl.rcParams["text.usetex"] = True
mpl.rcParams["font.family"]="serif"
mpl.rcParams["font.serif"]= "mathptmx"
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import stats

def LinearFit(x,y):
    slope,intercept,p_value,_,std_err = stats.linregress(x,y)
    return slope,intercept,std_err

def getPDF( field, nbins ):
        '''
        Outputs: 
                The normalized pdf of the input "field" over "nbins"
        '''
        field           = field.ravel()
        hist, bins  = np.histogram( field, bins=nbins )
        bins            = (bins[1:]+bins[:-1])*0.5
        pdf             = hist/(np.sum(hist)*np.mean(np.diff(bins)))

        return pdf, bins


data = pd.read_csv('data-send/shell_lam.csv')
print(data)
Re = data['Res']
Lambda = data['lam']
error = data['lam_err']
Slope_Shell,Intercept_Shell,Error_Shell = LinearFit(np.log(Re),np.log(Lambda))

print('SlopeShell=',Slope_Shell, 'Error Shell =',Error_Shell)

data = np.loadtxt('data-send/lambdaVsnu.txt')
Re_DNS = data[:,0]
Lambda_DNS = data[:,1]
error_DNS = data[:,2]
Slope_DNS,Intercept_DNS,Error_DNS = LinearFit(np.log(Re_DNS),np.log(Lambda_DNS))
data = np.loadtxt('data-send/gamma_avg_nu.txt')
gamma3 = data[:,0]
gamma3_std = data[:,3]
Slope_gamma3,Intercept_gamma3,Error_gamma3 = LinearFit(np.log(Re_DNS),np.log(-gamma3))
print('SlopeDNS=',Slope_DNS, 'Error DNS =',Error_DNS)
print('SlopeGamma3=',Slope_gamma3, 'Error Gamma3 =',Error_gamma3)

data = np.loadtxt('Post-Process-Data/LambdaMeanVsRe.txt')
LambdaMean = data[:,1]
ErrorLambdaMean = data[:,2]

slope, intercept, r_value, p_value, std_err = stats.linregress(np.log(Re),np.log(Lambda))
print('Slope Mean Lambda :',slope, 'Error:',std_err)


plt.rcParams.update({'font.size': 21})
fig, ax1 = plt.subplots(figsize=(8,7))
ms_size=10
m_ew=1.5 
alpha_=0.8
color_1='#014f86'  # Re DNS Color 
color_2='#bb4d00'  # Re Goy Color 
color_3 = '#2a6f97'
color_4 = '#2a6f97'

ax1.plot(Re_DNS, gamma3_std, 'd', ms=ms_size , fillstyle='none', color=color_1,label=r'$\gamma_3^{\rm std}$',mew=m_ew,alpha=1)
ax1.plot(Re_DNS,-gamma3,marker='^',ms=ms_size,fillstyle='full',ls='none',color=color_1,label=r'$-\langle \gamma_3 \rangle $',mew=m_ew,alpha=alpha_)
ax1.plot(Re_DNS,-gamma3[1]*(Re_DNS/Re_DNS[1])**0.5,ls='--',color=color_1,mew=m_ew,alpha=alpha_)

#ax1.plot(Re_DNS,Re_DNS**(Slope_gamma3)*np.exp(Intercept_gamma3),'--',color=color_3,mew=m_ew)

# ax1.errorbar(Re_DNS, Lambda_DNS, yerr=error_DNS,        fmt='s', ls='none', ms=ms_size, fillstyle='none', mew=m_ew,color=color_1,capsize=5,label=r'$\lambda_{\rm DNS}$')
# ax1.plot(Re_DNS,Re_DNS**(Slope_DNS)*np.exp(Intercept_DNS),'--',color=color_1,lw=2,alpha=alpha_)
ax1.errorbar(Re_DNS, LambdaMean, yerr=ErrorLambdaMean,  fmt='o', ls='none', mec='none',ms=ms_size,color=color_1,alpha=0.9,capsize=0,label=r'$ \lambda $')
ax1.errorbar(Re_DNS, Lambda_DNS,yerr=0, marker='s', ls='none', mec='none',ms=ms_size,color=color_1,alpha= 0.6,capsize=0,label=r'$\lambda^{\rm DNS}$')



ax1.set_xlabel(r'${\rm Re} (\rm{ DNS})$')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.spines['bottom'].set_color(color_1)
ax1.xaxis.label.set_color(color_1)
ax1.tick_params(axis='x', colors=color_1)
ax1.spines['bottom'].set_linewidth(2)
ax1.set_xlim([40,1e3])


ax2 = ax1.twiny()
ax2.errorbar(Re, 0.09*Lambda, yerr=error, fmt='^', ms=ms_size,  fillstyle='full', mec='none',mew=m_ew,capsize=0,color=color_2,label=r'$\lambda^{\rm GOY}$')
# ax2.plot(Re,0.18*Re**(Slope_Shell)*np.exp(Intercept_Shell),'--',color='k',lw=2,alpha=alpha_)
# print(np.exp(Intercept_Shell))

ax2.plot(Re,0.16*np.exp(Intercept_Shell)*Re**(0.59),'--',color='k',lw=2,alpha=alpha_)#,label=r'$Re^{0.58}$')
ax1.plot(Re_DNS,0.11*Re_DNS**(0.59),'--',color='k',lw=2,alpha=alpha_)#,label=r'$Re^{0.58}$')

print(Slope_DNS)
ax2.set_xlabel(r'${\rm Re} (\rm{GOY})$')
ax2.set_xscale('log')
ax2.spines['top'].set_color(color_2)
ax2.xaxis.label.set_color(color_2)
ax2.tick_params(axis='x', colors=color_2)
ax2.spines['top'].set_linewidth(2)
ax1.set_ylim([4*10**(-1),12])
lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()

# Combine and sort the labels
combined = sorted(zip(labels_1 + labels_2, lines_1 + lines_2), key=lambda x: x[0])
labels, lines = zip(*combined)
#plt.legend(ncol=2,fontsize=20)
ax1.legend(lines,labels, ncol=2, columnspacing=1 ,loc=[0.008,0.748], fontsize=18.5)

# Add text "Re" to the plot
ax1.text(300, 8, r'$\sim\mathrm{\sqrt{Re}}$', fontsize=18, color=color_1, rotation=0)
ax1.text(110,2.5, r'$\sim\mathrm{Re^{0.59}}$', fontsize=18, color='k', rotation=0)
plt.tight_layout()


# Inset 

nu = 0.016
'''
# S Lyapunov Exponent
data = np.loadtxt('Post-Process-Data/Local_Lyap_S_PDF_nu_%.3f.txt'%nu)
print(data.shape)
bins = data[0,:]
pdf = data[1,:]
binsp= bins[bins>0]
pdf_p = pdf[bins>0]
pdf_p = pdf_p
#pdf_p = pdf_p/np.sum(pdf_p)
bins_n = -bins[bins<0]
pdf_n = pdf[bins<0]
pdf_n = pdf_n
#pdf_n = pdf_n/np.sum(pdf_n)
ax_fontsize = 12
ax_inset = inset_axes(plt.gca(), width="50%", height="50%", loc='lower left', bbox_to_anchor=(0.40, 0.22, 0.38, 0.38), bbox_transform=plt.gcf().transFigure)
ax_inset.plot(binsp,pdf_p,'o' ,fillstyle='right',alpha=0.8,ms=4,label=r'$\lambda_{S}>0$')
ax_inset.plot(bins_n,pdf_n,'s',fillstyle='right',alpha=0.8,ms=4,label=r'$\lambda_{S}<0$')
# ax_inset.set_xscale('log')
# plt.legend(fontsize=ax_fontsize-5,loc=[0.10,0.15])
plt.legend(fontsize=ax_fontsize-4,loc='upper right')
ax_inset.set_yscale('log')
ax_inset.set_ylabel(r'$\mathcal{P}(\lambda_{S})$', fontsize=ax_fontsize)
ax_inset.set_xlabel(r'$\lambda_{S}$', fontsize=ax_fontsize)

# Set ticks
ax_inset.tick_params(axis='both', which='major', labelsize=ax_fontsize)
ax_inset.tick_params(axis='both', which='minor', labelsize=ax_fontsize)
# plt.yscale('log')
'''
# Visc Lyapunov Exponent in inset
ax_fontsize = 16
ax_inset = inset_axes(plt.gca(), width="50%", height="50%", loc='lower left', bbox_to_anchor=(0.58, 0.21, 0.61, 0.61), bbox_transform=plt.gcf().transFigure)
data = np.loadtxt('Post-Process-Data/Local_Lyap_Visc_PDF_nu_%.3f.txt'%nu)
print(data.shape)
bins = data[0,:]
pdf = data[1,:]
binsp= bins[bins>0]
pdf_p = pdf[bins>0]
pdf_p = pdf_p
#pdf_p = pdf_p/np.sum(pdf_p)
bins_n = -bins[bins<0]
pdf_n = pdf[bins<0]
pdf_n = pdf_n
#pdf_n = pdf_n/np.sum(pdf_n)

ax_inset.plot(binsp,pdf_p,'o' ,fillstyle='right',alpha=0.8,ms=4,label=r'$\lambda_{\eta}>0$')
ax_inset.plot(bins_n,pdf_n,'s',fillstyle='right',alpha=0.8,ms=4,label=r'$\lambda_{\eta}<0$')
ax_inset.plot(bins_n,bins_n**(-4),'--k',label=r'$\sim\lambda^{-4}$')
ax_inset.set_xscale('log')
# plt.legend(fontsize=ax_fontsize-8,loc=[0.25,0.58])
plt.legend(fontsize=ax_fontsize-1,loc='upper right')
ax_inset.set_yscale('log')
ax_inset.set_ylabel(r'$\mathcal{P}(\lambda_{\eta})$', fontsize=ax_fontsize)
ax_inset.set_xlabel(r'$\lambda_{\eta}$', fontsize=ax_fontsize)
ax_inset.set_xlim([20,10**3])
# Set ticks
ax_inset.tick_params(axis='both', which='major', labelsize=ax_fontsize)
ax_inset.tick_params(axis='both', which='minor', labelsize=ax_fontsize)





plt.savefig('LambdaVsReInset.jpg', dpi=200, bbox_inches='tight')
# plt.show()
