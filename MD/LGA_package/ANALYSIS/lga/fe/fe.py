import sys 
import numpy
import pylab 
import matplotlib.pyplot as plt 

dist=open('dist.xvg').readlines()
force=open('md_pull.xvg').readlines()

output=open('fe.out','w')

def read_file(file):
	coord_dict = dict()
	for line in file:
		if line.startswith('@') or line.startswith('#'):
			pass
		else:
			newline=line.split()
			x=float(newline[0])
			y=float(newline[1])
			coord_dict[x]=y
	return coord_dict 

dist_dict=read_file(dist)
force_dict=read_file(force)

for time in sorted(force_dict.keys()):
	if time in dist_dict:
		output.write(str(time) +'\t' + str(dist_dict[time]) +'\t' + str(force_dict[time])+ '\n')
		#print time, dist_dict[time], force_dict[time]

output.close()

extend=pylab.loadtxt('fe.out')
plt.plot(extend[:,1],extend[:,2])
plt.xlabel("Extension",fontsize="16")
plt.ylabel("Force", fontsize="16")
plt.show()
