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
			u_line = np.zeros((N+2))
			for i in range((N+2)):
				u_line[i] = float(infile.readline())
			u.append(u_line)
		x = np.linspace(0,u_line[-1], N+2)
		t = np.array(t)
		u = np.array(u)

		return x, t, u, dx, dt, Time, N

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
			u_line = np.zeros((N+2,N+2))
			for i in range(N+2):
				for j in range(N+2):

					u_line[i,j] = float(infile.readline())
					#print(u_line[i,j], i, j)

			ux.append(u_line[:,i])
			uy.append(u_line[i,:])


		x = np.linspace(0,u_line[-1], N+2)
		t = np.array(t)
		ux = np.array(ux)
		uy = np.array(uy)


		return x, t, ux, uy


def TwoDimPlot(x,y,z):
	X,Y = np.meshgrid(x,y)
	print(np.shape(z), np.shape(X))
	print(z[0])
	plt.pcolormesh(X,Y,z)
	plt.colorbar()
	plt.xlabel("x", fontsize=14)
	plt.ylabel("y", fontsize=14)
	plt.show()

	"""
	X,Y = np.meshgrid(x,y)
	print(np.shape(z[:,0:12]))
	for i in range(1,12):
		plt.pcolormesh(X,Y,z[:,12*i:12*(i+1)])
		plt.colorbar()
		plt.xlabel("x", fontsize=14)
		plt.ylabel("y", fontsize=14)
		plt.show()
	"""


def ThreeDimPlot(x,t,ux,uy):
	"""fig = plt.figure()
	print(np.shape(ux), np.shape(t))
	ax = fig.gca(projection="3d")
	ax.plot(ux[:,0,0], uy[:,0,0], t)
	plt.show()"""

	X, T = np.meshgrid(x, t)
	plt.pcolormesh(ux)
	plt.colorbar()
	plt.show()
	plt.plot(ux,uy)
	plt.show()



if __name__ == "__main__":
	filename = "main.txt"
	x, t, u, dx, dt, Time, N= read_dump(filename)
	x = np.linspace(0,1,144)
	y = np.linspace(0,1,26)
	print(np.shape(u))
	TwoDimPlot(x,y,u)
	#x, t, ux, uy = read_dump2(filename)
	#ThreeDimPlot(x,t,ux,uy)
