# super simple plotter..
import numpy as np
import matplotlib.pyplot as plt



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
			u_line = np.zeros(N+1)
			for i in range(N+1):

				u_line[i] = float(infile.readline())
			u.append(u_line)
		t = np.array(t)
		u = np.array(u)

		return t, u


# def TwoDimPlot(t,u):
#






if __name__ == "__main__":
	filename = "main.txt"
	t,u = read_dump(filename)
	plt.plot(u[5])
	plt.show()
