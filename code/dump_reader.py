import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

def read_dump2(filename):
	with open(filename, "r") as infile:
		N = int(infile.readline().split("=")[-1])
		Time = float(infile.readline().split("=")[-1])
		dx = float(infile.readline().split("=")[-1])
		dt = float(infile.readline().split("=")[-1])
		t = []
		ux = [] #Row
		uy = [] #Col
		for line in infile:
			t.append(float(line.split("=")[-1]))
			u_line = np.zeros((N+1,N+1))
			for i in range(N+1):
				for j in range(N+1):

					u_line[i,j] = float(infile.readline())
					#print(u_line[i,j], i, j)
			
			ux.append(u_line)
			#uy.append(u_line)

		x = np.linspace(0,u_line[-1], N+2)
		t = np.array(t)
		ux = np.array(ux[0])
		uy = np.array(uy[0])
		print(ux)

		return x, t, ux, uy


def TwoDimPlot(x,t,u):
	X,T = np.meshgrid(x,t)
	plt.pcolormesh(X,T,u)
	plt.colorbar()
	plt.xlabel("x", fontsize=14)
	plt.ylabel("time", fontsize=14)
	plt.show()


def ThreeDimPlot(x,t,ux,uy):
	fig = plt.figure()
	ax = fig.gca(projection="3d")

	ax.plot(ux)
	plt.show()




if __name__ == "__main__":
	filename = "main.txt"
	x, t, u = read_dump(filename)
	TwoDimPlot(x,t,u)
	#x, t, ux, uy = read_dump2(filename)
	#ThreeDimPlot(x,t,ux,uy)
