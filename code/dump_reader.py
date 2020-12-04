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




def read_dump(filename):
	with open(filename, "r") as infile:
		N = int(infile.readline().split("=")[-1])
		Time = float(infile.readline().split("=")[-1])
		dx = float(infile.readline().split("=")[-1])
		dt = float(infile.readline().split("=")[-1])
		t = []
		u = []
		for line in infile:
			t.append(float(line.split("=")[-1]))
			u_line = np.zeros(N+2)
			for i in range(N+2):
				u_line[i] = float(infile.readline())
			u.append(u_line)
		x = np.linspace(0,u_line[-1], N+2)
		t = np.array(t)
		u = np.array(u)

		return x, t, u


def TwoDimPlot(x,y,z):
	X,Y = np.meshgrid(x,y)
	plt.pcolormesh(X,Y,z)
	plt.colorbar()
	plt.xlabel("x", fontsize=14)
	plt.ylabel("y", fontsize=14)
	plt.show()










if __name__ == "__main__":
	filename = "main.txt"
	x, t, u = read_dump(filename)
	TwoDimPlot(x,t,u)
