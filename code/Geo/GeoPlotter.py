import sys
sys.path.append("..")
from dump_reader import *


if __name__ == "__main__":
    x, t, u, dx, dt, Time, N = read_dump(filename = "1DPlainDist_BRQ01.txt")
    OneDimMultiTime(x, t, u)
    #Animation_show(x, t, u)
