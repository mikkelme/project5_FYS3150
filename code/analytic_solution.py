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


def g(x,y,L=1):
    return np.sin(x)*np.sin(y)

def TwoDim_analytic(x, y, t_list, L=1):
    u = np.zeros((len(t_list), len(x), len(y)))
    for i in range(len(t_list)):
        t = t_list[i]
        for j in range(len(x)):
            u[i,j] = np.sin(np.pi*x[j]/L)*np.sin(np.pi*y/L)*np.exp(-(np.pi**2/L**2 + np.pi**2/L**2)*t)
    return u




    n_max = 10
    u = np.zeros((len(t_list), len(x), len(y)))
    for i in range(len(t_list)):
        t = t_list[i]
        for j in range(len(x)):
            sum = 0
            for n in range(1, n_max+1):
                npi = n*np.pi
                # B_x = L**2/np.pi**2*(n*np.sin(np.pi*L)*np.cos(npi) - L*np.cos(np.pi*L)*np.sin(npi))/(L**2 - n**2)
                #B_x = L/(L**2- npi**2)*(npi*np.sin(L)*np.cos(npi) - L*np.cos(L)*np.sin(npi))
                #B_x = L**2/(npi**2)*(np.sin(npi) - npi*np.cos(npi))
                # B_x = L*np.sin(npi)/(np.pi-np.pi*n**2)
                for m in range(1, n_max+1):
                    mpi = m*np.pi
                    #B_y = (m*np.sin(np.pi*L)*np.cos(mpi) - L*np.cos(np.pi*L)*np.sin(mpi))/(L**2 - m**2)
                    #B_y = L**2/(mpi**2)*(np.sin(mpi) - mpi*np.cos(mpi))
                    # B_y = (L-L*np.cos(npi))/npi
                    # B_y = L*np.sin(mpi)/(np.pi-np.pi*m**2)
                    # add = B_x*B_y*np.sin(npi*x[j]/L)*np.sin(npi*y/L)*np.exp(-(npi**2/L**2 + mpi**2/L**2)*t)
                    add = (1 + (-1)^(n+1))*(1-np.cos(mpi/2))/(n*m)*np.sin(npi/2*x[j])*np.sin(mpi/2*y)*np.exp(-np.pi**2*(n**2*m**2)*t/36)
                    sum += add
            u[i,j] = 200/np.pi**2*sum
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
