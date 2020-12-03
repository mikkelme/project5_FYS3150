import numpy as np
import matplotlib.pyplot as plt


def OneDim(x,t, L):
    n_max = 100
    v = 0
    for n in range(1, n_max):
        npi = n*np.pi
        add = 2*np.cos(npi)/(npi)*np.sin(npi*x/L)*np.exp(-npi**2*t/L**2)
        v += add
    u = v - f(x,L)
    return u

def f(x,L):
    return -x/L




x = np.linspace(0,1,100)
t = 0.00
L = 1
u = OneDim(x,t, L)
plt.plot(x,u)
plt.show()
