#! /usr/bin/env python

# FoldPlotLGA - version 5
# Author: Juami van Gils
# Last update: 14-03-2018
# ==========================
# Utilizes LGA and MD simulations to visualize the (un)folding curve of a given peptide/protein.
# Uses specific distance cutoff instead of absolute (average) GDTS to plot the unfolding landscape.

import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import os
import sys
import pylab
import ConfigParser
import scipy.optimize as optimization
from math import sqrt

#os.chdir(os.path.dirname(sys.argv[0]))
from savitzkygolay import savitzky_golay

#############################################################################
#FUNCTION: PARSE LGA
#INPUT: name of folder with frames, frame start, frame end, boolean: using thresholds to plot unfolding pathway?
#RETURNS: data matrix based on LGA files, list of residues
#############################################################################
def parse_lga(folder,frame_start, frame_end, binary_plot, column, treshold):
	#init empty data array	
	files = []
	#first iteration (file)
	firstFile = True
	residues = []
	#i = 1
	for frame in range(frame_start,frame_end+1):
		try:
			distances = []
			#when opening a new file, skip first line (header)
            #header = True
			file = '%sframe%s.pdb.frame0.pdb.res' % (folder, `frame`)
			with open(file) as input_data:
				for line in input_data:
					#start of block...
					line = line.strip()
					if line.startswith("LCS_GDT   NAME NUMBER"):
						break
				for line in input_data:
					#strip the line of linebreaks, split the line by space
					line = line.strip()
					dist = line.split()
					#check if line is empty, if TRUE then stop reading the file
					if line.startswith("LCS_AVERAGE"):
						break
					else:
						try:
							#folded
							if float(dist[column]) > treshold:
								distances.append(0)
							#unfolded binary
							elif binary_plot:
								distances.append(2)
							#unfolded colored
							else:
								distances.append(float(dist[column]))
							#list with residue numbers
							if firstFile:
								residues.append(dist[2])
						except:
							break
				#if not first file then stack array; otherwise create array
				#print distances
				if not firstFile:					
					data = np.row_stack((data, distances))
				#	i += 1
				#	print i
				else:
					data = np.array(distances)
					firstFile=False
				#if i == 50:
				#	break
		except:
			print "Skipping frame "+str(frame)
	return data, residues

#############################################################################
#FUNCTION: PARSE MD_PULL.XVG
#INPUT: path to md_pull.xvg
#RETURNS: dictionary with coordinates
#############################################################################
def parse_fe(file):
	coord_dict = dict()
	try:
		with open(file) as input_data:
			for line in input_data:
				if not line.startswith('@') and not line.startswith('#'):
					line = line.strip().split()
					coord_dict[round(float(line[0]), 2)] = float(line[1])
	except:
		print "ERROR: could not read file!"
	return coord_dict 

#############################################################################
#FUNCTION: PLOT FORCE-EXTENSION CURVE
#INPUT: dictionary force[time], pulling speed, initial length of the structure,
#		dictionary distance[time], boolean: analysis mode?, minimum force to plot, maximum force to plot
#RETURNS: boolean indicating whether or not the plot was created
#############################################################################
def plot_fe(force,pull_speed,init_len,time_dist, analysis_mode, min_force, max_force):
	fe_plotted = False
	try:
		print "INFO: trying to read md_pull.xvg..."
		force_dict = parse_fe(force)
		ax = plt.subplot(2,1,1)

		# get x and y coordinates from dictionaries
		time_force_dict = {}
		for key in time_dist:
			time_force_dict[time_dist[key]] = force_dict[key]
		x = sorted(time_force_dict.keys())
		y = [time_force_dict[k] for k in x]
		
		time_dist_rev = {}
		for key in time_dist:
			time_dist_rev[time_dist[key]] = key


		n_datapoints = len(x)
		
		# Smoothing
		smooth_window = n_datapoints / 160
		# Smoothing window must be an odd number
		if smooth_window % 2 == 0:
			smooth_window += 1
		
		# Step size for finding events
		step_size = n_datapoints / 40
		
		
		# Smooth curve using savitzky_golay method
		# window size given in input, polynomial order 3
		print 'Smoothing curve'
		y_array = np.asarray(y)
		yhat = savitzky_golay(y_array, smooth_window, 3) 

		# Find unfolding events
		event_indeces = []
		event_xcoords = []
		event_ycoords = []
		curve_starts = []

		print 'Searching for unfolding events'
		#int(round(1.5*step_size)
		for i in range(2*step_size, len(yhat)-step_size, step_size):
			cur_max = max(yhat[i:i+step_size])
			next_max = max(yhat[i+step_size:i+(2*step_size)])
			#if cur_max > next_max and not abs((cur_max - next_max)/cur_max) < 0.3:
			if cur_max > next_max and 7 < x[list(yhat).index(cur_max)] < 80:
				indx = list(yhat).index(cur_max)
				event_indeces.append(indx)
				event_xcoords.append(x[indx])
				event_ycoords.append(yhat[indx])

		# Add last event (only for full unfolding in simulations or
		# AFM experiments with ramp size around 80 nm)
		event_indeces.append(len(yhat)-1)
		event_xcoords.append(x[len(yhat)-1])
		event_ycoords.append(yhat[len(yhat)-1])

		# Start of fitting window is the minimum between two curves.
		# Skip the first step size of events for the first event, since
		# the start of the unfolding curve is usually quite noisy.
		for i in range(len(event_indeces)):
			if i == 0:
				curve_starts.append(2*step_size)
			else:
				curve_starts.append(list(yhat).index(min(yhat[event_indeces[i-1]:event_indeces[i]])))

		print 'Events found:', len(event_indeces)
		
		### WLC function
		def WLC(x, Lc):
			""" Function for WLC at T=300K, Lp=0.4. kb units converted to match x,Lp,Lc in nm and F in pN """
			Lp = 0.4
			kb = 1.38064852 * 10**(-2)				# SI: m^2 kg s^-2 K^-1
			return kb*300/Lp * (1/(4*(1 - x/Lc)**2) + x/Lc - 1/4)
			
		### Function to fit WLCs through curves
		def fit():
			""" Fit WLC through events and determine lcs """
			# New list to store fitted LC values and fitting errors
			LCs = []
			errs = []
			
			for i in range(len(event_indeces)):
				# Determine x window for WLC fitting and corresponding y-values
				window_start = curve_starts[i]
				window_end = event_indeces[i]
				xvals = x[window_start:window_end]
				yvals = yhat[window_start:window_end]
				
				# Estimate Lc as the x-coordinate of the unfolding event
				estimates = np.array(x[window_end])
				# Obtain Lp and Lc from curve fitting
				try:
					params, cov = optimization.curve_fit(WLC, xvals, yvals, estimates)
					Lc = params[0]
					LCs.append(Lc)
					err = sqrt(sum([(yvals[v] - WLC(xvals[v], Lc))**2 for v in range(len(xvals))]) / len(xvals))
					errs.append(err)
			
					# Add WLC to plot
					if i == len(event_indeces) - 1:
						xvals_plotting = x[:]
						xcurve = [time_dist_rev[v] for v in xvals_plotting]
					else:
						j = window_end
						x_new = x[j+1]
						f_new = WLC(x_new, Lc)
						while f_new < 800:
							x_new += 0.1
							f_new = WLC(x_new, Lc)
					
						for q in range(len(x)):
							if x[q] > x_new:
								j = q
								break		
						xvals_plotting = x[0:j]
						xcurve = [time_dist_rev[v] for v in xvals_plotting]
					plt.plot(xcurve, WLC(xvals_plotting, Lc), color='black', zorder=3)
				except:
					print 'Optimal fitting not found for range', min(xvals), max(xvals)	

			
			return event_indeces, LCs, event_xcoords, event_ycoords, errs
		
		# Fit WLC to curves
		event_indeces, LCs, event_xcoords, event_ycoords, errs = fit()
		
		xplot = sorted(time_dist_rev.values())
		yplot = [time_force_dict[time_dist[v]] for v in xplot]
		xevent = [time_dist_rev[v] for v in event_xcoords]
			
		plt.plot(xplot,yplot, zorder=0)
		plt.plot(xplot, yhat, color='red', zorder=1)
		plt.scatter(xevent, event_ycoords, color='black', zorder=2)
		plt.xlim(xmin=min(x), xmax=max(x))
		plt.ylim(ymin=min_force, ymax=max_force)
		plt.xlabel('End-to-end Distance (nm)', fontsize='16')
		plt.ylabel('Force (pN)', fontsize='16')
		ticks = ax.get_xticks()
		xaxis = ax.axes.get_xaxis()
		xaxis.set_ticks_position("top")
		
		
		#not in analysis mode, so plot altered axes
		plt.ylabel("Force (kJ/mol/nm)", fontsize="16")
		plt.xlim(xmin=min(force_dict.keys()), xmax=max(force_dict.keys()))
		plt.ylim(ymin=min_force, ymax=max_force)		
		ticks = ax.get_xticks()
		xaxis = ax.axes.get_xaxis()
		xaxis.set_ticks_position("top")
		#not in analysis mode, so plot altered axes
		plt.xlabel('End-to-end Distance (nm)', fontsize='16')
		ax.xaxis.set_label_position('top') 
		labels = []
		for label in ticks:
			if label == 0:
				labels.append(round(init_len,2))
			else:
				labels.append(round(time_dist[float(label)],2))
		ax.set_xticklabels(labels)
		
		plt.setp(ax.get_xticklabels(), fontsize=14)
		plt.setp(ax.get_yticklabels(), fontsize=14)
		fe_plotted = True
		print "INFO: plotting force-extension curve."
	except:
		print "WARNING: not plotting force-extension curve. Missing file(s)! Please check the 'fe' folder for 'md_pull.xvg' and 'dist.xvg'."
	
	return fe_plotted, event_indeces, LCs, event_xcoords, event_ycoords, errs

#############################################################################
#FUNCTION: PARSE DIST.XVG
#INPUT: path to file
#RETURNS: dictionary distance[time]
#############################################################################
def parse_dist_xvg(path):
	#dictionary with all [K][V] = [time][extension]
	time_dist = dict()
	with open(path) as dist_xvg:
		for line in dist_xvg:
			if not line.startswith('@') and not line.startswith('#'):
				#strip the line of linebreaks, split the line by spaces
				line = line.strip().split()			 
			  	#writing time/extension pairs to dictionary for easy access 
				time_dist[round(float(line[0]), 2)] = float(line[1])
	return time_dist

#############################################################################
#FUNCTION: PLOT BASED ON LGA DATA
#INPUT: lga data, location of output file, boolean: is the force-extension curve plotted,
#		frame start, frame end, pulling speed, frame step, pull start, initial length of the structure,
#		dictionary distance[time], boolean: analysis mode?, boolean: using thresholds to plot unfolding pathway?
#RETURNS: nothing
#############################################################################
def plot_lga(data,residues,outfile,fe_plotted,frame_start,frame_end,pull_speed,frame_step,fold_start,init_len,time_dist, analysis_mode, binary_plot, fig):
	print "INFO: plotting lga (un)folding map."
	#if statement: if fe plot has been made create a subplot (multipanel), otherwise just the one panel.
	v_min = 0
	v_max = max([max(x) for x in data])
	if fe_plotted:
		plot = plt.subplot(2,1,2)
		plt.xlim(xmax=frame_end/frame_step)
		if binary_plot:
			heatmap = plot.pcolor(data.transpose(), cmap=plt.cm.binary)
		else:
			heatmap = plot.pcolor(data.transpose(), cmap=plt.cm.afmhot, norm=colors.Normalize(vmin=0, vmax=25))
		
		if not analysis_mode:
			#convert labels: column number in matrix to extension
			xlabels = []
			for label in plot.get_xticks():
				if label == 0:
					xlabels.append(round(init_len,2))
				else:
					xlabels.append(round(time_dist[float(label*frame_step)],2))
			plot.set_xticklabels(xlabels)
					
			ylabels = []
			for label in plot.get_yticks():
				ylabels.append(int(label + fold_start))
			plot.set_yticklabels(ylabels)
			plt.xlabel("End-to-end Distance (nm)",fontsize="16")
			plt.ylabel("Residue", fontsize="16")
		else:
			plt.xlabel("Frame number",fontsize="16")
			plt.ylabel("Residue (relative)", fontsize="16")
			
			#end reset labels
	else:
		#reset plots
		plt.close()
		#plot unfolding curve
		if binary_plot:
			heatmap = plt.pcolor(data.transpose(), cmap=plt.cm.binary)
		else:
			heatmap = plt.pcolor(data.transpose(), cmap=plt.cm.afmhot)			
	
	plt.ylim(ymax=len(residues),ymin=1)
	plt.setp(plot.get_xticklabels(), fontsize=14)
	plt.setp(plot.get_yticklabels(), fontsize=14)
	cbaxes = fig.add_axes([0.92, 0.15, 0.02, 0.3]) 
	cbar = plt.colorbar(heatmap, cax = cbaxes)
	cbar.ax.tick_params(labelsize=14) 
	plt.savefig(outfile)

#############################################################################
#FUNCTION: CREATE PLOT (GENERAL)
#INPUT: ConfigParser object containing configuration file
#RETURNS: nothing
#############################################################################	
def create_plot(config):
	#get file related parameters
	path = config.get("file", "path")
	out = config.get("file", "out")
	stats = config.get("file", "stats")
	dist_xvg = config.get("file", "dist_xvg")
	pull_xvg = config.get("file", "pull_xvg")
	frame_start = int(config.get("simulation", "frame_start"))
	frame_end = int(config.get("simulation", "frame_end"))
	
	#get simulation parameters
	pull_speed = float(config.get("simulation", "pull_speed"))
	frame_step = int(config.get("simulation", "frame_step"))
	fold_start = int(config.get("simulation", "fold_start"))
	treshold = int(config.get('simulation', 'treshold'))
	cutoff = float(config.get('simulation', 'cutoff'))
	#create dictionary to convert the LGA distance cutoff to the corresponding column in the LGA output files
	LGA_dict = {}
	for i in list(np.arange(0.5, 10.5, 0.5)):
		LGA_dict[i] = int(2*i + 7)
	column = LGA_dict[cutoff]
	
	#get plot parameters
	analysis_mode = config.getboolean("plots","analysis")
	binary_plot = config.getboolean("plots","binary")
	min_force = int(config.get("plots", "min_force"))
	max_force = int(config.get("plots", "max_force"))

	#create dist[time] dictionary
	time_dist = parse_dist_xvg(dist_xvg)
	#print time_dist.keys()[:100]
	#exit()
	init_len = time_dist[0.00]
	
	fig = plt.figure(figsize=(8,6))
	data, residues = parse_lga(path,frame_start,frame_end, binary_plot, column, treshold)
	fe_plotted, event_indeces, LCs, event_xcoords, event_ycoords, errs = plot_fe(pull_xvg,pull_speed,init_len,time_dist, analysis_mode, min_force, max_force)
	plot_lga(data,residues,out,fe_plotted,frame_start,frame_end,pull_speed, frame_step,fold_start,init_len,time_dist, analysis_mode, binary_plot, fig)
	if binary_plot:
		print "INFO: Using thresholds to create unfolding pathway plot..."
	plt.subplots_adjust(hspace=0.05)
	plt.show()
	
	with open(stats, 'w') as f:
		f.write('Event_index\tBreak_extension\tLC\tBreak_force\tFitting_error\n')
		f.write('\n'.join(['%d\t%f\t%f\t%f\t%f' % (event_indeces[n], event_xcoords[n], LCs[n], event_ycoords[n], errs[n]) for n in range(len(LCs))]))

def main():
	config_file = sys.argv[1]
	config = ConfigParser.ConfigParser()
	config.read(config_file)
	create_plot(config)	

if __name__ == "__main__":
	sys.exit(main())
