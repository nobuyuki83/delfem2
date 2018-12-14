from OpenGL.GL import *

import sys
sys.path.append("../module_py")
import dfm2

import numpy, math

def main():

	A0 = numpy.zeros((128,128))
	for ix in range(A0.shape[0]):
		for iy in range(A0.shape[1]):		
			A0[iy,ix] = 10*math.sin(ix*0.1)*math.cos(iy*0.3)

	hm0 = dfm2.hight_map(A0)
	print(hm0.nDim,hm0.elemType)

	axis = dfm2.AxisXYZ()
	dfm2.winDraw3d([hm0,axis],(400,400))


if __name__ == "__main__":
  main()