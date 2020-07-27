from array_io import *
import numpy as np

fname = "HVF.txt"
fp = open(fname,"r")
fl = fp.readlines()
fp.close()

nred = 8
nv = 20
z = np.zeros(nred)

j = 0



for ired in range(nred):
	z[ired] = double(fl[j])
	j += 1