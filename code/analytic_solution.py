import numpy as np
import matplotlib.pyplot as plt

import sys
import seaborn as sns

plt.style.use("bmh")
sns.color_palette("hls", 1)

import matplotlib
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'



def OneDim_analytic(x, t_list, L=1):
    n_max = 100

    u = np.zeros((len(t_list), len(x)))
    for i in range(len(t_list)):
        v = 0
        t = t_list[i]
        for n in range(1, n_max):
            npi = n*np.pi
            add = 2*np.cos(npi)/(npi)*np.sin(npi*x/L)*np.exp(-npi**2*t/L**2)
            v += add
        u[i] = v - f(x,L)
    return u

def f(x,L):
    return -x/L




if __name__ == "__main__":
    x = np.linspace(0,1,100)
    t = [0, 1, 2]
    u = OneDim_analytic(x,t, L=1)
