#! /usr/bin/env python

"""
This script can be used to analyze all AFM curves in the input folder. It smooths the
data, corrects for drift, estimates the contact point of the tip and
gold surface, identifies unfolding events and fits the WLC model through
the events.
To run:
python AFM_analysis dir_in dir_out dir_out
"""

from sys import argv
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimization
import pylab
from math import ceil, sqrt
#from optparse import OptionParser, OptionGroup
curdir = os.getcwd()
os.chdir(os.path.dirname(argv[0]))
from savitzkygolay import savitzky_golay
os.chdir(curdir)

def y(x, a, b):
	""" Straight line with slope a and y-intersect b """
	return a*x + b
	
#def find_approach_retract_match(x, ar_coords):
#	""" for a given x value in the approach or retract
#	curve, find the index in the other curve that matches
#	that x value """
#	# Calculate differences
#	diff = [abs(x - n) for n in ar_coords]
#	m = min(diff)
#	idx = diff.index(m)
#	return idx

def smooth_data(approach, retract, outfile):
	""" Apply Savitzky-Golay smoothing on the curves """
	
	smooth_approach = list(savitzky_golay(np.asarray(approach[1]), 191, 3))
	smooth_retract = list(savitzky_golay(np.asarray(retract[1]), 191, 3))
	
	plt.figure()
	plt.plot(approach[0], approach[1], color='blue', zorder=0)
	plt.plot(retract[0], retract[1], color='red', zorder=1)
	plt.plot(approach[0], smooth_approach, color='green', zorder=2)
	plt.plot(retract[0], smooth_retract, color='green', zorder=2)
	plt.xlabel('Extension (nm)', fontsize='16')
	plt.ylabel('Force (pN)', fontsize='16')
	plt.title('Approach-Retract Curve with Smoothing', fontsize='18')
	plt.savefig(outfile)
	#plt.show()
	plt.close()
	
	
	return [approach[0], smooth_approach], [retract[0], smooth_retract]
	
def subtract_drift(approach, retract, drift_slope=None):
	""" Takes two lists containing lists with coordinates
	of the approach	and retract curves. returns the coordinates
	drift subtraction.
	"""
	
	# find the flat part of the approach curve
	start_fixed = 0
	step_size = len(approach[0])/15
	end = step_size
	
	# Find start of flat part
	estimates =  np.array([0, 0])
	params, cov = optimization.curve_fit(y, approach[0][start_fixed:end], approach[1][start_fixed:end], estimates)
	a_fit = params[0]
	b_fit = params[1]
	
	# Move to next segment until the slope is small
	while a_fit > 0.5:
		start_fixed += step_size
		end += step_size
		if end > len(approach[0])-1:
			start_fixed = 0
			break
		estimates =  np.array([a_fit, b_fit])
		params, cov = optimization.curve_fit(y, approach[0][start_fixed:end], approach[1][start_fixed:end], estimates)
		a_fit = params[0]
		b_fit = params[1]
	start = start_fixed
	
	# Move to next segment until slope increases again
	while a_fit < 3:
		start += step_size
		end += step_size
		if end > len(approach[0])-1:
			end = step_size
			break
		estimates =  np.array([a_fit, b_fit])
		params, cov = optimization.curve_fit(y, approach[0][start:end], approach[1][start:end], estimates)
		a_fit = params[0]
		b_fit = params[1]
	end -= step_size
	
	# Calculate drift. If drift_slope is given, use that as the slope of the fit
	if drift_slope == None:
		estimates =  np.array([0, b_fit])
		params, cov = optimization.curve_fit(y, approach[0][start_fixed:end], approach[1][start_fixed:end], estimates)	
		drift_slope = params[0]
		b_fit = params[1]
		print 'drift', drift_slope, b_fit
	else:
		estimates =  np.array([b_fit])
		params, cov = optimization.curve_fit(lambda x, b: y(x, drift_slope, b), approach[0][start_fixed:end], approach[1][start_fixed:end], estimates)	
		b_fit = params[0]
		print 'drift', drift_slope, b_fit

	# Calculate new force values from approach and retract curves	
	y_approach = [approach[1][n] - y(approach[0][n], drift_slope, b_fit) for n in range(len(approach[0]))]
	y_retract = [retract[1][n] - y(retract[0][n], drift_slope, b_fit) for n in range(len(retract[0]))]

	"""
	# Plot corrected curves
	print drift_slope, b_fit
	yvals = [y(n, drift_slope, b_fit) for n in approach[0][start_fixed:end]]
	plt.figure()
	plt.plot(approach[0], approach[1], color='blue', zorder=0)
	plt.plot(retract[0], retract[1], color='red', zorder=1)
	plt.plot(approach[0][start_fixed:end], yvals, color='green', zorder=2)
	plt.show()
	"""
	return [approach[0], y_approach], [retract[0], y_retract], drift_slope

def flip_axes(approach, retract):
	""" Flip axes such that the extension increasees with x and the
	event forces become positive """
	
	ax_flip = [-n for n in approach[0]]
	ay_flip = [-n for n in approach[1]]
	rx_flip = [-n for n in retract[0]]
	ry_flip = [-n for n in retract[1]]
	"""
	plt.figure()
	plt.plot(ax_flip, ay_flip, color='blue', zorder=0)
	plt.plot(rx_flip, ry_flip, color='red', zorder=1)
	plt.show()
	"""
	return [ax_flip, ay_flip], [rx_flip, ry_flip]

def k_xcorr(approach, retract, k=None):
	""" Correct the x coordinates and if necessary determine spring constant"""
	
	# Take points up to a force smaller than -80 pN
	for n in retract[1]:
		if not n < -80:
			end = retract[1].index(n)
			break
	print 'end', end
	if k == None:
		estimates =  np.array([0, 0])
		params, cov = optimization.curve_fit(y, retract[0][:end], retract[1][:end], estimates)	
		a_fit = float(params[0])
		b_fit = float(params[1])
	else:
		estimates =  np.array([0])
		params, cov = optimization.curve_fit(lambda x, b: y(x, k, b), retract[0][:end], retract[1][:end], estimates)	
		a_fit = k
		b_fit = float(params[0])
	
	# Find the contact point of the curve (y=0 of the fitted line)
	# and set this x value to zero
	x_zero = -b_fit/a_fit
	approach[0] = [v - x_zero for v in approach[0]]
	retract[0] = [v - x_zero for v in retract[0]]
	
	
	"""
	yvals = [y(n, a_fit, b_fit) for n in retract[0][:end]]
	plt.figure()
	plt.plot(approach[0], approach[1], color='blue', zorder=0)
	plt.plot(retract[0], retract[1], color='red', zorder=1)
	plt.plot(retract[0][:end], yvals, color='black', zorder=2)
	plt.show()
	"""
	return approach, retract, a_fit

def bending_correction(approach, retract, k, f_name, d_out):
	""" Correct curve for tip bending """
	
	## Fit line through the first part of the curve where the
	## tip pushes into the glass curve. Look at deviation from
	## this line to determine where to start searching for events
	#for n in retract[1]:
	#	if not n < -80:
	#		end = retract[1].index(n)
	#		break
	
	#estimates =  np.array([0, 0])
	#params, cov = optimization.curve_fit(y, retract[0][:end], retract[1][:end], estimates)	
	#a_fit = float(params[0])
	#b_fit = float(params[1])
	"""
	yvals = [y(n, a_fit, b_fit) for n in retract[0][:end]]
	plt.figure()
	plt.plot(approach[0], approach[1], color='blue', zorder=0)
	plt.plot(retract[0], retract[1], color='red', zorder=1)
	plt.plot(retract[0][:end], yvals, color='black', zorder=2)
	plt.show()
	"""
	# Find the contact point of the curve (y=0 of the fitted line)
	# and set this x value to zero
	#x_zero = -b_fit/a_fit
	#approach[0] = [v - x_zero for v in approach[0]]
	#retract[0] = [v - x_zero for v in retract[0]]
	
	# spring constant is the slope of the fitted line
	#k = a_fit
	
	x_approach = []
	for i in range(len(approach[0])):
		x_approach.append(approach[0][i] - approach[1][i]/k)
	
	x_retract = []
	for i in range(len(retract[0])):
		x_retract.append(retract[0][i] - retract[1][i]/k)
	"""	
	plt.figure()
	plt.plot(x_approach, approach[1], color='blue', zorder=0)
	plt.plot(x_retract[end:], retract[1][end:], color='red', zorder=1)
	plt.plot(x_retract[:end], retract[1][:end], color='black', zorder=2)
	plt.show()
	"""
		
	outfile = os.path.join(d_out, 'approach_' + f_name)
	with open(outfile, 'w') as f:
		f.write('\n'.join(['%f\t%f' % (x_approach[i], approach[1][i]) for i in range(len(approach[1]))]))
	
	outfile = os.path.join(d_out, 'retract_' + f_name)
	with open(outfile, 'w') as f:
		f.write('\n'.join(['%f\t%f' % (x_retract[i], retract[1][i]) for i in range(len(retract[1]))]))
	
	return [x_approach, approach[1]], [x_retract, retract[1]]

def find_events(approach, retract):
	""" Find unfolding events in the retract curve """
	# Minimal drop in force to call an event
	event_size = 20
	# Minimum extension value
	alpha = 5
	# Maximum force to start counting events (to get rid of adhesion)
	beta = 2000
	# Segmentation window
	step_size = 250
	# Number of points to fit force drop
	drop_size = 200
	# Start at end if no good start point can be found (== skip curve)
	start = len(retract[0]) - 1
	
	# Find the first point where the extension is larger than
	# the offset
	print "Don't take into account points at x lower than", alpha
	for i in range(len(retract[0])):
		if retract[0][i] > alpha:
			start = i
			break
	
	# Find first point after start with retract close enough to zero
	# force to avoid fitting adhesion
	print 'Find first point with force lower than', beta, 'to avoid fitting adhesion'
	for i in range(start, len(retract[1])):
		if retract[1][i] < beta:
			start = i
			break
	
	# Find unfolding events
	print 'Find unfolding events'
	event_indeces = []
	
	for i in range(start, len(retract[0]) - step_size, step_size):
		cur_max = max(retract[1][i:i+step_size])
		next_max = max(retract[1][i+step_size:i+(2*step_size)])
		#print 'max', cur_max, next_max
		# Event if segment drops with at least event size
		if cur_max - next_max > 2:
			indx = i + list(retract[1][i:i+step_size]).index(cur_max)
			
			# Fit line after first drop_size values after the maximum value of the segment.
			# An event is called if the slope is very negative or if a force drop larger than
			# event_size is present in this region
			estimates =  np.array([0, 0])
			params, cov = optimization.curve_fit(y, retract[0][indx:indx+drop_size], retract[1][indx:indx+drop_size], estimates)	
			a_fit2 = float(params[0])
			b_fit2 = float(params[1])
			print 'a_fit', a_fit2
			if (a_fit2 < -15 or cur_max - min(retract[1][indx:indx+drop_size]) > event_size):# and cur_max > 14:
				event_indeces.append(indx)
	
	# If the maximum in the last segment of the retract curve has a force larger than 
	# event_size call that maximum as well
	last_max = max(retract[1][-step_size:])
	if last_max > event_size:
		idx = retract[1].index(last_max)
		if not idx in event_indeces and idx > start:
			event_indeces.append(idx)
	
	return event_indeces, start
	
def fit_WLC(approach, retract, start, event_indeces, f_name):
	""" Fit WLC through unfolding events """
	# Output file names
	outfile1 = f_name + '_fit.png'
	outfile2 = f_name + '_zoom_fit.png' 
	# Max number of points for WLC fitting
	fitting_size = 2000
	# Number of points to fit force drop
	drop_size = 200
	
	# Make figure
	plt.figure()
	plt.plot(approach[0], approach[1], color='blue', zorder=0)
	plt.plot(retract[0], retract[1], color='red', zorder=1)
	
	# Initiate lists to store relevant variables
	curve_starts = []
	event_xcoords = []
	event_ycoords = []
	force_drops = []
	
	for i in range(len(event_indeces)):
		indx = event_indeces[i]
		# max fitting window of 500 datapoints
		min_index = event_indeces[i]-fitting_size
		if i == 0:
			curve_starts.append(max(start, min_index))
			#plt.plot(retract[0][max(start, min_index):event_indeces[i]], retract[1][max(start, min_index):event_indeces[i]], color='green', zorder=2)
			#plt.plot(retract[0][indx:indx+drop_size], retract[1][indx:indx+drop_size], color='yellow', zorder=2)
		else:
			# Find the minimum between two events
			between_events = retract[1][event_indeces[i-1]:event_indeces[i]]
			min_middle = event_indeces[i-1] + between_events.index(min(between_events))
			curve_starts.append(max(min_middle, min_index))
			#plt.plot(retract[0][max(min_middle, min_index):event_indeces[i]], retract[1][max(min_middle, min_index):event_indeces[i]], color='green', zorder=2)
			#plt.plot(retract[0][indx:indx+drop_size], retract[1][indx:indx+drop_size], color='yellow', zorder=2)
	
	# Remove first point if start and end of fitting window are the same
	for i in range(len(event_indeces)):
		if event_indeces[i] == curve_starts[i]:
			event_indeces.remove(event_indeces[i])
			curve_starts.remove(curve_starts[i])
			break
	
	for indx in event_indeces:
		event_xcoords.append(retract[0][indx])
		event_ycoords.append(retract[1][indx])
		force_drops.append(retract[1][indx] - min(retract[1][indx:indx+drop_size]))
	
	# remove events that are too close together (and therefore probably belong to the same event)
	i = 0
	while i < len(event_indeces) - 1:
		for i in range(len(event_indeces)):
			print i
			xval = event_xcoords[i]
			if len([x for x in event_xcoords if xval-4 < x < xval+4]) > 1:
				print 'Double event found'
				event_indeces.remove(event_indeces[i+1])
				curve_starts.remove(curve_starts[i+1])
				event_xcoords.remove(event_xcoords[i+1])
				event_ycoords.remove(event_ycoords[i+1])
				force_drops.remove(force_drops[i+1])
				break
		print event_indeces
	
	"""
	plt.figure()
	plt.plot(approach[0], approach[1], color='blue', zorder=0)
	plt.plot(retract[0][:-step_size], retract[1][:-step_size], color='red', zorder=1)
	plt.plot(retract[0][-step_size:], retract[1][-step_size:], color='black', zorder=1)
	plt.scatter(event_xcoords, event_ycoords, color='black', zorder=2)
	plt.show()
	"""
	def WLC(x, Lc):
		""" Function for WLC at T=300K, Lp=0.4. kb units converted to match x,Lp,Lc in nm and F in pN """
		Lp = 0.4
		kb = 1.38064852 * 10**(-2)		# SI: m^2 kg s^-2 K^-1
		T = 300
		return kb*T/Lp * (1/(4*(1 - x/Lc)**2) + x/Lc - 1/4)
		
		
	# List to store fitted LC values
	LCs = []
	errs = []
	plt.scatter(event_xcoords, event_ycoords, color='black', zorder=2)

	print 'Fit WLC to the curve'
	for i in range(len(event_indeces)):
		# Determine x window for WLC fitting and corresponding y-values
		window_start = curve_starts[i]
		window_end = event_indeces[i]
		xvals = retract[0][window_start:window_end]
		yvals = retract[1][window_start:window_end]
		print window_start, window_end
		# Estimate Lc as slightly larger than the x-coordinate of the unfolding event
		estimates = np.array(retract[0][window_end]+1)
		weights = np.linspace(2, 1, len(xvals))
		#print estimates
		# Obtain Lc from curve fitting
		#try:
		params, cov = optimization.curve_fit(WLC, xvals, yvals, estimates, sigma=weights, absolute_sigma=True)
		Lc = params[0]
		# sum of squares error divided by the number of datapoints used for the fitting
		#err = sum([(retract[1][v] - WLC(retract[0][v], Lc))**2 for v in range(event_indeces[i]-300:event_indeces[i])])
		err = sqrt(sum([(yvals[v] - WLC(xvals[v], Lc))**2 for v in range(len(xvals))]) / len(xvals))
		errs.append(err)
		LCs.append(Lc)
		
		# Add WLC to plot
		xvals_plotting = retract[0][:]
			
		j = window_end
		x_new = retract[0][min(j+1, len(retract[0])-1)]
		f_new = WLC(x_new, Lc)
		while f_new < 700 and x_new < 150:
			x_new += 0.5
			f_new = WLC(x_new, Lc)
			#print 'f', x_new, f_new, Lc
	
		for q in range(len(retract[0])):
			if retract[0][q] > x_new:
				xvals_plotting = retract[0][:q]
				break
			elif q == len(retract[0]) - 1:
				xvals_plotting += range(int(ceil(retract[0][-1])), int(ceil(x_new)), 1)
		#print 'q', q, retract[0][q]
		#print xvals_plotting[-1]
		plt.plot(xvals_plotting, WLC(xvals_plotting, Lc), color='black', zorder=3)
	# Set plotting parameters
	plt.xlim(xmin=-10, xmax=max(retract[0])+10)
	plt.ylim(ymin=-800, ymax=max(retract[1])+50)
	plt.xlabel('Extension (nm)', fontsize='16')
	plt.ylabel('Force (pN)', fontsize='16')
	plt.title('Misfolded', fontsize='18')
	plt.savefig(outfile1)
	#plt.show()
	plt.close()
	
	#plot zoomed in
	plt.figure()
	plt.plot(approach[0], approach[1], color='blue', zorder=0)
	plt.plot(retract[0], retract[1], color='red', zorder=1)
	
	for i in range(len(event_indeces)):
		indx = event_indeces[i]
		#plt.plot(retract[0][curve_starts[i]:event_indeces[i]], retract[1][curve_starts[i]:event_indeces[i]], color='green', zorder=2)
		#plt.plot(retract[0][indx:indx+drop_size], retract[1][indx:indx+drop_size], color='yellow', zorder=2)
		
		Lc = LCs[i]
		j = event_indeces[i]
		xvals_plotting = retract[0][:]
			
		x_new = retract[0][min(j+1, len(retract[0])-1)]
		f_new = WLC(x_new, Lc)
		while f_new < 700 and x_new < 150:
			x_new += 0.1
			f_new = WLC(x_new, Lc)
			#print 'f', x_new, f_new, Lc
	
		for q in range(len(retract[0])):
			if retract[0][q] > x_new:
				xvals_plotting = retract[0][:q]
				break
			elif q == len(retract[0]) - 1:
				xvals_plotting += range(int(ceil(retract[0][-1])), int(ceil(x_new)), 1)
		plt.plot(xvals_plotting, WLC(xvals_plotting, Lc), color='black', zorder=3)
	
	plt.scatter(event_xcoords, event_ycoords, color='black', zorder=2)
	
	plt.xlim(xmin=-10, xmax=max(retract[0])+10)
	plt.ylim(ymin=-50, ymax=max(retract[1])+50)
	plt.xlabel('Extension (nm)', fontsize='16')
	plt.ylabel('Force (pN)', fontsize='16')
	plt.title('Misfolded', fontsize='18')
	plt.savefig(outfile2)
	plt.close()
	
	return LCs, event_xcoords, event_ycoords, force_drops, errs
	
	
def main():
	d_in = argv[1]
	d_out = argv[2]
	corr_out = argv[3]
	
	if not os.path.isdir(d_out):
		os.mkdir(d_out)
	
	if not os.path.isdir(corr_out):
		os.mkdir(corr_out)
	else:
		print 'Directory', corr_out, 'already exists'
		#exit(-1)
	
	unique_curves = [f for f in os.listdir(d_in) if 'retract' in f]
	
	spring_constants = []
	# Estimate k from all retraction curves
	for retract_file in unique_curves:
				
		r_x = []
		r_y = []
		
		# Read approach curve coordinates
		with open(os.path.join(d_in, retract_file), 'r') as f:
			for line in f:
				line = line.strip().split('\t')
				r_x.append(float(line[0]))
				r_y.append(float(line[1]))
			
		r_x = [-v for v in r_x]
		r_y = [-v for v in r_y]
	
		estimates =  np.array([0, 0])
		params, cov = optimization.curve_fit(y, r_x[:500], r_y[:500], estimates)	
		k_fit = float(params[0])
		spring_constants.append(k_fit)
		
	k = float(np.median(spring_constants))
	
	found_coords = []
	for f_name in unique_curves:
		
		split_name = f_name.split('_')
		print split_name
		f_x = split_name[1]
		f_y = split_name[2]
		
		if not (f_x, f_y) in found_coords:
			retract_curves = sorted([f for f in unique_curves if f_x in f and f_y in f])
			approach_curves = sorted(['approach_' + f.split('retract_')[1] for f in retract_curves])
			print '# curves:', len(approach_curves), len(retract_curves)
			
			for i in range(len(approach_curves)):
				outname = approach_curves[i].split('approach_')[1]
				
				# Read in coordinate files
				approach = os.path.join(d_in, approach_curves[i])
				retract = os.path.join(d_in, retract_curves[i])
				
				
				a_x = []
				a_y = []
				r_x = []
				r_y = []
				
				# Read approach curve coordinates
				with open(approach, 'r') as f:
					for line in f:
						line = line.strip().split('\t')
						a_x.append(float(line[0]))
						a_y.append(float(line[1]))
					
				# Read retract curve coordinates
				with open(retract, 'r') as f:
					for line in f:
						line = line.strip().split('\t')
						r_x.append(float(line[0]))
						r_y.append(float(line[1]))
				
				a = [a_x, a_y]
				r = [r_x, r_y]
				
				# Smoothing
				print 'Smoothing curve'
				outfile = os.path.join(d_out, outname + '_smooth.png')
				a_smooth, r_smooth = smooth_data(a, r, outfile)
				try:
					# Subtract drift
					print 'Subtracting drift'
					if i == 0:
						a_drift, r_drift, drift_slope = subtract_drift(a_smooth, r_smooth)
					else:
						a_drift, r_drift, drift_slope = subtract_drift(a_smooth, r_smooth, drift_slope)
					print min(r_drift[1]), max(r_drift[1]), drift_slope
					# Flip axes
					print 'Flip axes to make extension and forces positive'
					a_flip, r_flip = flip_axes(a_drift, r_drift)

					# Determine k and correct for x = 0
					print 'Determining spring constant and correct for x=0'
					if not k == None:
						a_xcorr, r_xcorr, k = k_xcorr(a_flip, r_flip, k)
					elif i == 0:
						a_xcorr, r_xcorr, k = k_xcorr(a_flip, r_flip)
					else:
						a_xcorr, r_xcorr, k = k_xcorr(a_flip, r_flip, k)
						
					print 'k', k
					
					# Find unfolding events
					event_indeces, start = find_events(a_xcorr, r_xcorr)
					
					# Correct for tip bending and write corrected coordinates into output file
					a_bend, r_bend = bending_correction(a_xcorr, r_xcorr, k, outname, corr_out)
					
					# Fit WLC through events
					LCs, x_breaks, forces, force_drops, errs = fit_WLC(a_bend, r_bend, start, event_indeces, os.path.join(d_out, outname))
					
					# Write LCs into output file
					outfile = os.path.join(d_out, outname + '_LCs.txt')
					with open(outfile, 'w') as f:
						f.write('Break_extension\tLC\tBreak_force\tForce_drop\tFitting_error\tMax_extension\n')
						f.write('\n'.join(['%f\t%f\t%f\t%f\t%f\t%f' % (x_breaks[n], LCs[n], forces[n], force_drops[n], errs[n], max(r_bend[0])) for n in range(len(LCs))]))
				except:
					print 'Failed for curve', outname
			found_coords.append((f_x, f_y))
			
if __name__ == '__main__':
	main()
