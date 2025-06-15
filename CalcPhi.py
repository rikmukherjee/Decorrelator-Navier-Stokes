# %%
import numpy as np
import h5py
import math
import array
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import ticker
from matplotlib import rcParams
from scipy.io import FortranFile as ff
from scipy.fftpack import fftn, ifftn




# %%
# Define the shape of the array
N = 64  
# Read the Fortran file
t = 2
folder = '%d/R2/vel/'%N

def read_velocities(t, folder):
    file_path = folder + 'u%d.in' % t
    with ff(file_path, 'r') as f:
        u = f.read_reals(dtype=np.float64).reshape((N, N, N))
    file_path = folder + 'v%d.in' % t
    with ff(file_path, 'r') as f:
        v = f.read_reals(dtype=np.float64).reshape((N, N, N))
    file_path = folder + 'w%d.in' % t
    with ff(file_path, 'r') as f:
        w = f.read_reals(dtype=np.float64).reshape((N, N, N))
    return u, v, w

u, v, w = read_velocities(t, folder)
print(u.shape, v.shape, w.shape)

# %%
def check_divergence(u, v, w):
    u_hat = fftn(u)
    v_hat = fftn(v)
    w_hat = fftn(w)

    # Create the wavenumber arrays (kx, ky, kz)
    nx, ny, nz = u.shape
    kx = np.fft.fftfreq(nx) * nx
    ky = np.fft.fftfreq(ny) * ny
    kz = np.fft.fftfreq(nz) * nz

    # Create the wavenumber meshgrid
    Kx, Ky, Kz = np.meshgrid(kx, ky, kz, indexing='ij')

    # Compute the strain rate tensor components in Fourier space
    Sxx = np.real(ifftn(1j * Kx * u_hat))  
    Syy = np.real(ifftn(1j * Ky * v_hat))
    Szz = np.real(ifftn(1j * Kz * w_hat))
    print('Divergence =', (Sxx + Syy + Szz).max())
    del Sxx, Syy, Szz




def Calc_Sij(u,v,w):
    u_hat = fftn(u)
    v_hat = fftn(v)
    w_hat = fftn(w)

    # Create the wavenumber arrays (kx, ky, kz)
    nx, ny, nz = u.shape
    kx = np.fft.fftfreq(nx) * nx
    ky = np.fft.fftfreq(ny) * ny
    kz = np.fft.fftfreq(nz) * nz

    # Create the wavenumber meshgrid
    Kx, Ky, Kz = np.meshgrid(kx, ky, kz, indexing='ij')
    Sxx = np.real(ifftn(1j * Kx * u_hat))  
    Syy = np.real(ifftn(1j * Ky * v_hat))
    Szz = np.real(ifftn(1j * Kz * w_hat))
    Sxy = np.real(ifftn(1j * 0.5 * (Ky * u_hat + Kx * v_hat)))
    Sxz = np.real(ifftn(1j * 0.5 * (Kz * u_hat + Kx * w_hat)))
    Syz = np.real(ifftn(1j * 0.5 * (Kz * v_hat + Ky * w_hat)))
    del u_hat, v_hat, w_hat

    return Sxx, Syy, Szz, Sxy, Sxz, Syz

# %%
Phi_t = []

for t in range(30):
    folder =  '%d/R2/vel/'%N
    u, v, w = read_velocities(t, folder)

    folder= '%d/R3/vel/'%N
    u_pert, v_pert, w_pert = read_velocities(t, folder)

    print('KE=', 0.5*(u**2+v**2+w**2).mean())
    phi =  0.5*((u_pert-u)**2 + (v_pert-v)**2 + (w_pert-w)**2)
    Phi_t.append(np.sum(phi/N**3))

    plt.savefig('Plots/phi%d.jpg' % t,bbox_inches='tight')
    plt.clf()




