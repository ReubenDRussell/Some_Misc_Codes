# Written by Reuben.
####################################
print("Clacks Overhead: GNU Terry Pratchett");
####################################
import math
import matplotlib.pyplot as plt
import numpy as np
import os
####################################


# ---------------------------------
# Cartesian to polar and back again
# ---------------------------------
def cart2pol(x, y):
	R = np.sqrt(x**2 + y**2);
	phi = np.arctan2(y, x);
	return(R, phi)

def pol2cart(R, phi):
	x = R * np.cos(phi);
	y = R * np.sin(phi);
	return(x, y)
# ---------------------------------

# ---------------------------------
# Spherical Coordinate Conversion
# ---------------------------------
def cart2sph(x,y,z):
	R = np.sqrt(x**2 + y**2 + z**2);
	theta = np.arccos(z/R);
	if(x>0):
		phi = np.arctan(y/x);
	if(x<0 and y>=0):
		phi = np.arctan(y/x) + np.pi;
	if(x<0 and y<0):
		phi = np.arctan(y/x) - np.pi;
	if(x==0 and y>0):
		phi = 0.5*np.pi;
	if(x==0 and y<0):
		phi = -0.5*np.pi;
	if(x==0 and y==0):
		return(print("Undefined theta due to x and y intersection locations both being 0. Breaking loop ... try again?"));
	return(R, theta, phi)

def sph2cart(R, theta, phi):
	x = R*np.sin(theta)*np.cos(phi);
	y = R*np.sin(theta)*np.sin(phi);
	z = R*np.cos(theta);
	return(x,y,z)

# ---------------------------------



def circle_intersection(x_orig, y_orig, radius_orig, x_offs, y_offs, radius_offs):
	dist_1_2 = np.sqrt((x_offs - x_orig)**2 + (y_offs - y_orig)**2);	# distance between circle 1 (original) and 2 (offset)
	#
	if(dist_1_2 > radius_orig + radius_offs):					# if distance greater than combind radii --> leave
		return []
	if(dist_1_2 < abs(radius_orig - radius_offs)):					# if one circle inside the other --> leave
		return []
	if(dist_1_2 == 0 and radius_orig == radius_offs):				# if circle 1 is circle 2 --> leave
		return []
	if(dist_1_2 == radius_orig + radius_offs or dist_1_2 == abs(radius_orig - radius_offs)):	# if circle 1 and circle 2 tangential --> do below (may need original to be 0,0 centre?) (copied from old version as I cant be bothered)
		P_x = x_orig + (x_offs - x_orig) * radius_orig / dist_1_2;			# x position of overlap
		P_y = y_orig + (y_offs - y_orig) * radius_orig / dist_1_2;			# y position of overlap
		return [(P_x, P_y)]
	# The equation for finding the x and y overlaps of a circle can be broken into two non-x-y terms to help simplify the code
	term1 = (radius_orig**2 - radius_offs**2)/(2*(dist_1_2**2));
	term2 = 0.5*np.sqrt( 2*((radius_orig**2 + radius_offs**2)/(dist_1_2**2)) - (((radius_orig**2 - radius_offs**2)**2)/(dist_1_2**4)) -1 );
	X_intersection1 = (0.5*(x_orig + x_offs)) + (term1*(x_offs - x_orig)) + (term2*(y_offs - y_orig));
	Y_intersection1 = (0.5*(y_orig + y_offs)) + (term1*(y_offs - y_orig)) + (term2*(x_orig - x_offs));
	X_intersection2 = (0.5*(x_orig + x_offs)) + (term1*(x_offs - x_orig)) - (term2*(y_offs - y_orig));
	Y_intersection2 = (0.5*(y_orig + y_offs)) + (term1*(y_offs - y_orig)) - (term2*(x_orig - x_offs));
	
	Intersection1 = (X_intersection1, Y_intersection1);
	Intersection2 = (X_intersection2, Y_intersection2);
	return[Intersection1, Intersection2]

	# FULL Equation of intersection points of two circles (the equation allowes for two different radii and central points)
	# X_intersection1 = 0.5*(x_orig + x_offs) + (( (radius_orig**2 - radius_offs**2)/(2*(dist_1_2**2))) * (x_offs - x_orig) ) + 0.5*(y_offs - y_orig)*np.sqrt( 2*((radius_orig**2 + radius_offs**2)/(dist_1_2**2)) - (((radius_orig**2 - radius_offs**2)**2)/(dist_1_2**4)) -1 )
	# Y_intersection1 = 0.5*(y_orig + y_offs) + (( (radius_orig**2 - radius_offs**2)/(2*(dist_1_2**2))) * (y_offs - y_orig) ) + 0.5*(x_orig - x_offs)*np.sqrt( 2*((radius_orig**2 + radius_offs**2)/(dist_1_2**2)) - (((radius_orig**2 - radius_offs**2)**2)/(dist_1_2**4)) -1 )
	# X_intersection2 = 0.5*(x_orig + x_offs) + (( (radius_orig**2 - radius_offs**2)/(2*(dist_1_2**2))) * (x_offs - x_orig) ) - 0.5*(y_offs - y_orig)*np.sqrt( 2*((radius_orig**2 + radius_offs**2)/(dist_1_2**2)) - (((radius_orig**2 - radius_offs**2)**2)/(dist_1_2**4)) -1 )
	# Y_intersection2 = 0.5*(y_orig + y_offs) + (( (radius_orig**2 - radius_offs**2)/(2*(dist_1_2**2))) * (y_offs - y_orig) ) - 0.5*(x_orig - x_offs)*np.sqrt( 2*((radius_orig**2 + radius_offs**2)/(dist_1_2**2)) - (((radius_orig**2 - radius_offs**2)**2)/(dist_1_2**4)) -1 )


def print_and_plot_intersections(radii, offset_x, offset_y, z_distance):
	if os.path.exists("Ring_Intersection_points.txt"):	# An important check on if an output file for this already exists. If not, it will continue.
		print("Rename or delete Ring_Intersection_points.txt before trying again");
		return()
	else:
		print("Ring_Intersection_points.txt doesn't exist, continuing operations ...");
		all_intersections = [];
		### Loop that goes through and utilises the intersection finding code for all rings which intersect
		for rad1 in radii:
			for rad2 in radii:
				intersections = circle_intersection(original_x, original_y, rad1, offset_x, offset_y, rad2);	# 0,0 is centre coord for unperturbed circle
				all_intersections.extend(intersections);
				# Printing intersection points to terminal
				if intersections:
					if len(intersections) == 1:  # tangential	(just about impossible unless a VERY big beam offset)
						print(f"Radii {rad1}, {rad2} -> Intersection: {intersections[0]}");
					else:
						x_test1 , y_test1 = intersections[0]	# These two added to extract the x and y individually from intersections arrays
						x_test2 , y_test2 = intersections[1]	#	(this way we can use them to see the polar coordinate output)
						R_from_targ1 , Theta1 , Phi1 = cart2sph(x_test1, y_test1, z_distance);
						R_from_targ2 , Theta2 , Phi2 = cart2sph(x_test2, y_test2, z_distance);
						with open("Ring_Intersection_points.txt", "a") as output_file:	# to to output it to the file aswell as the terminal
							print(f"{rad1} , {rad2}\t\t{R_from_targ1} , {Theta1} , {Phi1}\t\t{R_from_targ2} , {Theta2} , {Phi2}", file = output_file);
							# print(f"{rad1} , {rad2}\t\t{R_from_targ1} , {np.degrees(Theta1)} , {np.degrees(Phi1)}\t\t{R_from_targ2} , {np.degrees(Theta2)} , {np.degrees(Phi2)}", file = output_file);
						# print(f"{rad1} , {rad2}\t || \t{R_from_targ1} , {Theta1} , {Phi1}\t || \t{R_from_targ2} , {Theta2} , {Phi2}"); # \t means insert tab
						print(f"{rad1} , {rad2}\t || \t{R_from_targ1} , {np.degrees(Theta1)} , {np.degrees(Phi1)}\t || \t{R_from_targ2} , {np.degrees(Theta2)} , {np.degrees(Phi2)}");
						# print(f"Radii {rad1}, {rad2} -> Spherical Intersections: {R_from_targ1}, {Theta1}, {Phi1} || {R_from_targ2}, {Theta2}, {Phi2}")
						# print(f"Radii {rad1}, {rad2} -> Cartesian Intersections: {intersections[0]}, {intersections[1]}")
						# print(f"Radii {rad1}, {rad2} -> Polar Intersections: {cart2pol(x_test1 , y_test1)}, {cart2pol(x_test2 , y_test2)}")

	####
	# Plotting
	###
	fig, ax = plt.subplots(figsize=(8, 8));	# How big do you want the figure to be?
	ax.set_xlim(-Figure_size, Figure_size);			# X and Y limits used to be +/- 40 for both. Now they are adaptive
	ax.set_ylim(-Figure_size, Figure_size);			# 
	ax.axhline(0, color='grey', lw=0.5);
	ax.axvline(0, color='grey', lw=0.5);
	### Plot the circles
	for r in radii:
		circle = plt.Circle((original_x, original_y), r, fill=False, color='red', linestyle='dashed', label='Circle 0,0');
		ax.add_artist(circle);
		circle_off = plt.Circle((offset_x, offset_y), r, fill=False, color='blue', linestyle='dashed', label='Offset Circle');
		ax.add_artist(circle_off);
	### Plot the intersection points
	ix, iy = zip(*all_intersections);
	ax.scatter(ix, iy, color='green', s=50, label='Intersection Points');
	### Plot the group start and end points
	orig_start_circ = plt.Circle((original_x, original_y), Group_Start_radius, fill=False, color='red', linestyle='solid', label='Circle 0,0', lw=0.5);	# solid line at start of ring group	original
	ax.add_patch(orig_start_circ);
	orig_end_circ = plt.Circle((original_x, original_y), Group_End_radius, fill=False, color='red', linestyle='solid', lw=0.5);				# ------"------ end of ring group	original
	ax.add_patch(orig_end_circ);
	offs_start_circ = plt.Circle((offset_x, offset_y), Group_Start_radius, fill=False, color='blue', linestyle='solid', label='Offset Circle', lw=0.5);	# solid line at start of ring group	offset
	ax.add_patch(offs_start_circ);
	offs_end_circ = plt.Circle((offset_x, offset_y), Group_End_radius, fill=False, color='blue', linestyle='solid', lw=0.5);				# ------"------ end of ring group	offset
	ax.add_patch(offs_end_circ);
	### Axes and final touches
	# print(f"TEST; x = {ix} , y = {iy}");
	plt.xlabel('X (mm)');
	plt.ylabel('Y (mm)');
	plt.title('Intersection Points of an S3 Detector Rings - offset vs original');
	plt.legend();
	plt.grid(True);
	plt.show();


###################################################################################
# Input required information here (Given data)
###################################################################################
# num_circles = 25;	# needs to be total number of rings +1 because it will only go to the beginning of the last ring (just in case)
# initial_radius = 11;	# initial radius (hole in middle of micron S3) is 11mm
##
num_circles = 8;	# for testing above and below the high and low rings
initial_radius = 9;	#	--------------""-------------
##
increment = 1;		# width of a ring (1mm)
z_distance = 33.0;	# the distance of the targets from the detectors
original_x = 0;
original_y = 0;
offset_x = -1.35;
offset_y = +0.65;
# offset_x = -1;
# offset_y = -1;
# Exp S1576 approx Beam offset:	(check signs are correct and make logical sense)
#	x= +1.35
#	y= -0.65
Group_Start = 1;	# start ring of group (remember index starts from 1)
Group_End = 4;		# end of ring group
# Group_Start_radius = Group_Start + initial_radius - increment;	# from start of first ring in the group
Group_Start_radius = 11; # if doing a test for below start ring "rings"
# Group_End_radius = Group_End + initial_radius;		# to the end of the last ring in the group
Group_End_radius = Group_End + Group_Start_radius;	# needs to be group made from start ring therefore Group_End is the group size?	# NOTE: the end radius is at the end of ring 24?
Figure_size = num_circles + initial_radius + 5;
###################################################################################
radii = [initial_radius + i * increment for i in range(0, num_circles, 1)];	# array of radii in the range defined by number of circles and initial radius i starts at 0 and steps by 1
print_and_plot_intersections(radii, offset_x, offset_y, z_distance);