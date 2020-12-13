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
		x = np.linspace(0,dx*(N+1), N+2)
		t = np.array(t)
		u = np.array(u)

		return x, t, u, dx, dt, Time, N

def Animation_show(x, t, u):
	import matplotlib.animation as animation
	fig = plt.figure()
	ax = plt.axes(xlim=(x[0], x[-1]), ylim = [np.min(u),np.max(u)])
	numerical, = ax.plot([], [], lw=2, label = "numerical")

	# initialization function: plot the background of each frame
	def init():
	    numerical.set_data([], [])
	    return numerical,

	# animation function.  This is called sequentially
	def animate(i):
	    i = int(i)
	    numerical.set_data(x,u[i])
	    ax.set_title(f"t = {t[i]} ({t[i]/t[-1]*100:.2f})%")
	    #ax.set_ylim(0, np.max(u_ana))
	    return numerical,

	# call the animator.
	time_frames = np.linspace(0,len(t)-1,len(t))
	anim = animation.FuncAnimation(fig, animate, init_func=init,
	                               frames=time_frames, interval=1, blit=False)

	plt.legend(bbox_to_anchor=(0.4, 0.80))
	plt.show()


def read_dump2(filename):
	with open(filename, "r") as infile:
		N = int(infile.readline().split("=")[-1])
		Time = float(infile.readline().split("=")[-1])
		dx = float(infile.readline().split("=")[-1])
		dt = float(infile.readline().split("=")[-1])
		t = []
		ux = [] #Row
		uy = [] #Col
		u = []
		for line in infile:
			t.append(float(line.split("=")[-1]))
			u_line = np.zeros((N+2,N+2))
			for i in range(N+2):
				for j in range(N+2):

					u_line[i,j] = float(infile.readline())
			u.append(u_line)


		xy = np.linspace(0, dx*(N+1), N+2)
		t = np.array(t)
		u = np.array(u)
		return xy, t, u, dx, dt, Time, N


def TwoDimSubplots(xy, t, u, dx, dt, Time, N):
	X,Y = np.meshgrid(xy,xy)
	colormap = "plasma"
	vmin = 0; vmax = 1
	print(dx, dt)

	fig = plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
	t_idx = np.linspace(0,len(t)-1, 4).astype(int)
	t_idx = [0,30,35,37]


	for i in range(4):
	 	plt.subplot(2,2,i+1)
	 	plt.title(f"t = {t[t_idx[i]]:.2f}")
	 	mesh = plt.pcolormesh(X,Y,u[t_idx[i]], vmin = vmin, vmax = vmax, cmap = colormap)

	plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
	fig.subplots_adjust(right = 0.8)
	cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
	fig.colorbar(mesh, cax = cbar_ax, label = "$u$")
	plt.show()




def TwoDimPlot(x,y,z):
	X,Y = np.meshgrid(x,y)

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

	x, t, u, dx, dt, Time, N = read_dump(filename = "1DPlainDist2.txt")
	Animation_show(x, t, u)


	# x, t, u, dx, dt, Time, N = read_dump(filename = "1DExplicit0.1.txt")
	# x = np.linspace(0,1,144)
	# y = np.linspace(0,1,26)

	# xy, t, u, dx, dt, Time, N = read_dump2(filename = "2DExplicit0.05.txt")
	# TwoDimSubplots(xy, t, u, dx, dt, Time, N)


	#TwoDimPlot(xy,xy,u[0])
	# ThreeDimPlot(x,t,ux,uy)
