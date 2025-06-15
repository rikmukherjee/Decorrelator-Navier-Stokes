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

N = 256
ReVals = np.array([692.78,358.25,178.53,98.72,49.58])

markers = ['o','s','^','D','v']
colors = ['#124559','#bb3e03','#f48c06','#faa307','#ffba08']
colors = ["#03045e","#023e8a","#0077b6","#0096c7","#00b4d8","#48cae4","#90e0ef","#ade8f4","#caf0f8"]
# colors = ["#dad7cd","#a3b18a","#588157","#3a5a40","#344e41"]
# colors = ["#6a040f","#9d0208","#d00000","#dc2f02","#e85d04","#f48c06","#faa307","#ffba08"] [2::] # Red
colors = ["#005f73","#0a9396","#ee9b00","#ca6702","#bb3e03","#ae2012"]
plt.figure(figsize=[8,6])
mpl.rcParams.update({'font.size': 26})
for nu_idx in range(1,6):
    data =np.loadtxt('data-send/Phi_nu_%d'%nu_idx)
    t   = data[:,0]
    Phi = data[:,1]
    # Logarithmically placed t index from 0 to len(t)
    log_idx = np.logspace(0,np.log10(len(t)-1),40,base=10,dtype=int)
    log_idx = np.unique(log_idx)


    #plt.plot(t[log_idx], Phi[log_idx],ls='-',lw=2, marker = markers[nu_idx - 1], color=colors[nu_idx-1],ms=10, fillstyle='right',label=r'$%.1f$'%ReVals[nu_idx-1])
    plt.plot(t[1::2], Phi[1::2],ls='-',lw=2, marker = markers[nu_idx - 1], color=colors[nu_idx-1],
             ms=10, fillstyle='right',label=r'$%.1f$'%ReVals[nu_idx-1])
    # plt.plot(t[1::2], Phi[1::2],ls='-',lw=2, marker = markers[nu_idx - 1],mew=2,
    #          color=colors[nu_idx-1],ms=10, fillstyle='none',label=r'$%.1f$'%ReVals[nu_idx-1],alpha=0.8)

plt.axhline(y=Phi[-1],ls='--',color='black',lw=2)
# plt.text(3,Phi[-1]*0.05,r'$2E_0$',fontsize=20)
plt.yscale('log')
plt.xlabel(r'$t$')
plt.ylabel(r'$\Phi(t)$')
plt.xlim(0,55)
plt.legend(title=r'$Re=$', title_fontsize=18,loc=[0.8,0.5],fontsize=14)
plt.tight_layout()


# Create Inset 
ax_inset = inset_axes(plt.gca(), width="50%", height="40%", loc='lower left', bbox_to_anchor=(0.4, 0.180, 0.72, 0.72), bbox_transform=plt.gcf().transFigure)
t=1
phi = np.load('256/phi%d.npy'%t)
ax_inset.imshow(phi[:,:,N//2],cmap='hot')
ax_inset.set_title(r'$t = %.1f$'%((t-1)/2),fontsize=20)
plt.axis('off')

ax_inset = inset_axes(plt.gca(), width="50%", height="40%", loc='lower right', bbox_to_anchor=(0.32, 0.180, 0.72, 0.72), bbox_transform=plt.gcf().transFigure)
t=13
phi = np.load('256/phi%d.npy'%t)
ax_inset.imshow(phi[:,:,N//2],cmap='hot')
ax_inset.set_title(r'$t = %.1f$'%((t-1)/2),fontsize=18)
plt.axis('off')

plt.savefig('PhiVs_t.jpg',dpi=300,bbox_inches='tight')
# plt.show()
















