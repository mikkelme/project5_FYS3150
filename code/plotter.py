# super simple plotter..
import numpy as np
import matplotlib.pyplot as plt


file = open("main.txt")
lines = file.readlines()
N = len(lines)

u = []
for i in range(N):
	u.append(float(lines[i]))

plt.plot(np.linspace(0,1,N),u)
plt.show()
