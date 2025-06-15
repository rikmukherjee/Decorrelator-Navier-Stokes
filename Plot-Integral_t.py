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
data = np.loadtxt('data-send/lambdaVsnu.txt')
lambdaVals = data[:,1]

LambdaMean = []
Error = []

for nu_idx in range(1,6):
    Lambda = lambdaVals[nu_idx-1]

    data = np.loadtxt('Post-Process-Data/Integral_t_nu_%d'%nu_idx)
    t   = data[0:-1,0]
    Int_Sij = data[0:-1,1]
    Int_Visc = data[0:-1,2]
    # Phi = (data[:,3]**2)/2
    data =np.loadtxt('data-send/Phi_nu_%d'%nu_idx)
    Phi = data[:,1]
    tIni = np.where(Phi>10**(-14))[0][0]
    tFin = np.where(Phi>10**(-5))[0][0]
    PhiDotByPhi = (-Int_Sij+ Int_Visc)/Phi
    Error.append(np.std(PhiDotByPhi[tIni:tFin]))
    LambdaMean.append(np.mean(PhiDotByPhi[tIni:tFin]))



    plt.plot(t, Int_Sij/Phi,ls='-',lw=2, label=r'$\int_{\Omega} S_{ij} d\Omega$',color='blue')
    plt.plot(t, Int_Visc/Phi,ls='-',lw=2, label=r'$\int_{\Omega} \nu \nabla u d\Omega$',color='red')
    plt.plot(t, PhiDotByPhi ,ls='-',lw=2, label=r'$\Phi(t)$',color='green')
    plt.axhline(y=Lambda,ls='--',color='black',lw=2)
    plt.axvline(x=t[tIni],ls='--',color='black',lw=2)
    plt.axvline(x=t[tFin],ls='--',color='black',lw=2)
    # plt.ylim(-4,4)
    plt.show()

LambdaMean = np.array(LambdaMean)
Error = np.array(Error)
savedata = np.column_stack((ReVals,LambdaMean,Error))
np.savetxt('Post-Process-Data/LambdaMeanVsRe.txt',savedata)














