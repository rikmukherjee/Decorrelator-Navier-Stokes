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
mpl.rcParams.update({'font.size': 25})
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



# plt.figure(figsize=(8,6))

# data = np.loadtxt('data-send/lambdaVsnu.txt')
# Re_DNS = data[:,0]
# data = np.loadtxt('data-send/gamma_avg_nu.txt')
# gamma3 = data[:,0]
# gamma3_err = data[:,3]
# slope,intercept,err = LinearFit(np.log(Re_DNS),np.log(gamma3_err))
# print('SlopeGamma3=',slope, 'Error Gamma3 =',err)

# plt.plot(Re_DNS,gamma3_err,'o',fillstyle='none',mew=2.2,ms=9)
# plt.plot(Re_DNS,gamma3_err[1]*(Re_DNS/Re_DNS[1])**(0.59),'--k')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel(r'$Re$')
# plt.ylabel(r'$Std(\gamma_3)$')
# plt.tight_layout()
# plt.show()

data = np.loadtxt('data-send/lambdaVsnu.txt')
Re_DNS = data[:,0]
markers = ['o','s','^','D','v']
for nu_idx in range(1,6):
        print(nu_idx)
        data = np.loadtxt('Post-Process-Data/Dist-Gamma3-nu-%d.txt'%nu_idx)
        bins= data[:,0]
        pdf = data[:,1]
        plt.plot(bins,pdf,label=r'${\rm Re=} = %.1f$'%Re_DNS[nu_idx-1],marker=markers[nu_idx-1],mfc='none',ms=5,fillstyle='right',mew=1.2)
plt.xlim(-50,0)
plt.ylabel(r'$\mathcal{P}(\gamma_3)$')
plt.xlabel(r'$\gamma_3$')
plt.legend(fontsize=14,loc='upper left')
plt.ylim(5*10**(-5),1)
plt.yscale('log')
plt.tight_layout()
plt.savefig('Distribution-Gamma3.png',dpi=300)
plt.show()