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

def f(x,L):
    return -x/L

def OneDim_analytic(x, t_list, L=1):
    n_max = 100
    u = np.zeros((len(t_list), len(x)))
    for i in range(len(t_list)):
        v = 0
        t = t_list[i]
        for n in range(1, n_max+1):
            npi = n*np.pi
            add = 2*np.cos(npi)/(npi)*np.sin(npi*x/L)*np.exp(-npi**2*t/L**2)
            v += add
        u[i] = v - f(x,L)
    return u


def TwoDim_analytic(x, y, t_list, L=1):
    u = np.zeros((len(t_list), len(x), len(y)))
    for i in range(len(t_list)):
        t = t_list[i]
        for j in range(len(x)):
            u[i,j] = np.sin(np.pi*x[j]/L)*np.sin(np.pi*y/L)*np.exp(-(np.pi**2/L**2 + np.pi**2/L**2)*t)
    return u








if __name__ == "__main__":
    x = np.linspace(0,1,100)
    y = np.linspace(0,1,100)
    t = [0, 0.2, 1]
    # u = OneDim_analytic(x,t, L=1)
    u = TwoDim_analytic(x,y,t)
    from dump_reader import *
    TwoDimPlot(x,y,u[2])























#
