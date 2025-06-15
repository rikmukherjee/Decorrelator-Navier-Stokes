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

data = pd.read_csv('data-send/shell_lam.csv')
print(data)
Re = data['Res']
Lambda = data['lam']
error = data['lam_err']
Slope_Shell,Intercept_Shell,Error_Shell = LinearFit(np.log(Re),np.log(Lambda))



data = np.loadtxt('data-send/lambdaVsnu.txt')
Re_DNS = data[:,0]
Lambda_DNS = data[:,1]
error_DNS = data[:,2]
Slope_DNS,Intercept_DNS,Error_DNS = LinearFit(np.log(Re_DNS),np.log(Lambda_DNS))

data = np.loadtxt('data-send/gamma_avg_nu.txt')
gamma3 = data[:,0]
Slope_gamma3,Intercept_gamma3,Error_gamma3 = LinearFit(np.log(Re_DNS),np.log(-gamma3))


data = np.loadtxt('Post-Process-Data/LambdaMeanVsRe.txt')
LambdaMean = data[:,1]
ErrorLambdaMean = data[:,2]


plt.rcParams.update({'font.size': 26})
fig, ax1 = plt.subplots(figsize=(8,7))
ms_size=10
m_ew=1.5
alpha_=0.8
color_1='#014f86'  # Re DNS Color 
color_2='#bb4d00'  # Re Goy Color 
color_3 = '#2a6f97'
color_4 = '#2a6f97'

ax1.plot(Re_DNS,- gamma3, 'd', ms=ms_size, fillstyle='full', color=color_1,label=r'$\gamma_3$',mew=m_ew,alpha=alpha_)
ax1.plot(Re_DNS,Re_DNS**(Slope_gamma3)*np.exp(Intercept_gamma3),'--',color=color_3,mew=m_ew)

# ax1.errorbar(Re_DNS, Lambda_DNS, yerr=error_DNS,        fmt='s', ls='none', ms=ms_size, fillstyle='none', mew=m_ew,color=color_1,capsize=5,label=r'$\lambda_{\rm DNS}$')
# ax1.plot(Re_DNS,Re_DNS**(Slope_DNS)*np.exp(Intercept_DNS),'--',color=color_1,lw=2,alpha=alpha_)
ax1.errorbar(Re_DNS, LambdaMean, yerr=ErrorLambdaMean,  fmt='o', ls='none', mec='none',ms=ms_size,color=color_1,alpha=0.9,capsize=0,label=r'$\langle \lambda \rangle$')
ax1.errorbar(Re_DNS, Lambda_DNS,yerr=0, marker='s', ls='none', mec='none',ms=ms_size,color=color_1,alpha= 0.6,capsize=0,label=r'$\lambda_{\rm DNS}$')



ax1.set_xlabel(r'${\rm Re}^{\rm{ DNS }}$')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.spines['bottom'].set_color(color_1)
ax1.xaxis.label.set_color(color_1)
ax1.tick_params(axis='x', colors=color_1)
ax1.spines['bottom'].set_linewidth(2)
ax1.set_xlim([40,1e3])


ax2 = ax1.twiny()
ax2.errorbar(Re, 0.2*Lambda, yerr=error, fmt='o', ms=ms_size,  fillstyle='full', mec='none',mew=m_ew,capsize=0,color=color_2,label=r'$\lambda_{\rm GOY}$')
# ax2.plot(Re,0.18*Re**(Slope_Shell)*np.exp(Intercept_Shell),'--',color='k',lw=2,alpha=alpha_)
ax2.plot(Re,0.22*Re**(0.58)*np.exp(Intercept_Shell),'--',color='k',lw=2,alpha=alpha_)#,label=r'$Re^{0.58}$')
print(Slope_DNS)
ax2.set_xlabel(r'${\rm Re}^{\rm{GOY}}$')
ax2.set_xscale('log')
ax2.spines['top'].set_color(color_2)
ax2.xaxis.label.set_color(color_2)
ax2.tick_params(axis='x', colors=color_2)
ax2.spines['top'].set_linewidth(2)

lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()

# Combine and sort the labels
combined = sorted(zip(labels_1 + labels_2, lines_1 + lines_2), key=lambda x: x[0])
labels, lines = zip(*combined)

ax1.legend(lines, labels, loc='upper left', fontsize=20)
plt.tight_layout()

plt.savefig('LambdaVsRe.jpg', dpi=200, bbox_inches='tight')
plt.show()
