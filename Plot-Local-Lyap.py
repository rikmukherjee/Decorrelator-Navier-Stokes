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


nuVals = 2**np.arange(1,6)*10**(-3)
ReVals = np.array([692.78,358.25,178.53,98.72,49.58])
plt.figure(figsize=[8,6])
mpl.rcParams.update({'font.size': 20})

nu = 0.016

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

plt.plot(binsp,pdf_p,'o' ,fillstyle='right',alpha=0.9,ms=8,label=r'$\lambda_{S}>0$')
plt.plot(bins_n,pdf_n,'s',fillstyle='right',alpha=0.7,ms=8,label=r'$\lambda_{S}<0$')
plt.plot(bins_n,10*(bins_n)**(-4),'--k',label=r'$\sim \lambda_{S}^{-4}$')
plt.ylim(10**(-10),10**(-1))
plt.xlim(10**(1),10**(3))
plt.legend(loc='lower left')
plt.xlabel(r'$\lambda_{S}$')
plt.ylabel(r'$\mathcal{P}(\lambda_{S})$')
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()

# Visc Lyapunov Exponent in inset
ax_fontsize = 16
ax_inset = inset_axes(plt.gca(), width="50%", height="50%", loc='lower left', bbox_to_anchor=(0.55, 0.55, 0.65, 0.65), bbox_transform=plt.gcf().transFigure)
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
# ax_inset.set_xscale('log')
plt.legend(fontsize=ax_fontsize-3,loc=[0.5,0.65])
ax_inset.set_yscale('log')
ax_inset.set_ylabel(r'$\mathcal{P}(\lambda_{\eta})$', fontsize=ax_fontsize)
ax_inset.set_xlabel(r'$\lambda_{\eta}$', fontsize=ax_fontsize)

# Set ticks
ax_inset.tick_params(axis='both', which='major', labelsize=ax_fontsize)
ax_inset.tick_params(axis='both', which='minor', labelsize=ax_fontsize)














plt.savefig('local_lyap_dist.png')
plt.show()



