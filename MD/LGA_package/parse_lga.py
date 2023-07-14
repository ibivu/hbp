#! /usr/bin/env python

import sys
import ntpath
import matplotlib
import matplotlib.pyplot as plt

distance = list()#init distance list
LGA_file= sys.argv[1]#define input file: first argument 
filename = ntpath.basename(LGA_file) #get name of file

print "Creating plot..."

with open(LGA_file) as input_data:
	# Skips text before the beginning of the interesting block:
	for line in input_data:
        	if line.strip() == '#      Molecule1      Molecule2  DISTANCE    Mis    MC     All    Dist_max   GDC_mc  GDC_all':  #start of block... NB: assuming that this is indeed the structure of the file c.q. header and this doesn't change as the input files change.
            		break
    	for line in input_data:  # end of block..
		dist = line.strip().split() #strip the line of linebreaks, split the line by spaces
		if len(dist) == 0: #check if line is empty, if TRUE then stop reading the file
			break	
        	distance.append(dist[5])  #add distance per residue to list

#write distance/residues to file:
with open("ANALYSIS/output/"+filename+".out","w") as outfile:
	res = 1 #residue number
	outfile.write("residue \t distance\n")
	for dist in distance:#write all distances to file
		outfile.write(str(res)+"\t"+dist+"\n")
		res += 1

#plotting...
plt.plot(distance)
plt.ylabel('Distance',fontsize=16)
plt.xlabel('Residue',fontsize=16)
plt.title('Residue-distance plot',fontsize=18)
plt.xlim(xmax=len(distance)-1) #set limit to x-axis, only need to plot as many points as there are residues

plt.savefig("ANALYSIS/output/"+filename+".out.png")

print "Done!"
