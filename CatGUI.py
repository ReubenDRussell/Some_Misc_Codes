import tkinter as tk;			# for GUI
import numpy as np;			# for maths
import math;				# for maths
# import IsotopeMassInfo as IsotopeInfo;			# custom script to read masses.txt
# import physics_constants as physc;	# custom script (specific physics constants and unit converters)	only needed for u_to_MeV?
# import physics_formulas as physf;	# custom script	(other useful physics equations)			UNUSED?
import pandas as pd;			# No fucking clue - black and white fluffy thing (not cuddly) -- array thing?

import datetime;			# for the date and time in the output txt file

import os;				# needed for the reading of the masses20.txt

import matplotlib as plt;		# for the plots of the main data
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk);	# specific extra bits
from matplotlib.figure import Figure;

## needed. Otherwise it cuts off the number I can show
pd.set_option('display.max_rows', None);
pd.set_option('display.max_columns', None);
pd.set_option('display.width',1000)

u_to_MeV       = 931.49432         # MeV / c^2	## Only thing taken from physics_constants.py ==> no need to import anymore


###############
# Definitions & Equations:
###############


#==============
# the below chunk is what was once in IsotopeMassInfo.py
#==============

list_element_sym =	["n","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K",
			"Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb",
			"Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs",
			"Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta",
			"W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa",
			"U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt",
			"Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"];	# updated - not bothering to make it use the txt table for symbols and Z yet

def get_Z( sym ):
	#Returns Z for given chemical symbol.
	listsym = list(sym);
	listsym[0] = listsym[0].upper();
	sym = ''.join(listsym);
	return list_element_sym.index(sym);

def get_sym( Z ):
	#Returns chemical symbol for given Z.
	tmpZ = int(Z);
	if tmpZ < 0 or tmpZ > len(list_element_sym):
		raise ValueError( 'Z={} is an invalid element number'.format(Z) );
	return list_element_sym[tmpZ];

def decompose_atomic_mass( sam,sam_err ):
	"""Decomposes a string of format 1.234(56)# into two parts: 1.234 and
	56. Returns float(1.234), float(56) and True if there are no # in
	the argument.  Use with caution! Requires really this format as
	input!
	"""
	ipo = sam.find('#');
	ipc = sam_err.find('#');
	iextrapolated = sam.find('#');
	iextrapolated_err = sam_err.find('#');
	if "#" in sam:
		sam1 = float(sam[0:ipo])/(1000000);	# because it is in micro-u and want in u (atomic mass units)
	else:
		sam1 = float(sam)/(1000000);

	if "#" in sam_err:
		sam_err1 = float(sam_err[0:ipc])/(1000000);
	else:
		sam_err1 = float(sam_err)/(1000000);

	return (float(sam1),	# float the number between column 0 and column ipo of input
		float(sam_err1),
		iextrapolated==-1, iextrapolated_err==-1); # adds a flag to these that they are extrapolated experimental values?

# Initialization of the IsotopeInfo module
def _initialize_module():
	# Getting the path where this module IsotopeInfo.py is located
	__location__ = os.path.realpath( os.path.join(os.getcwd(), os.path.dirname(__file__)));
	# print('Initializing module IsotopeInfo')
	with open( os.path.join(__location__, 'masses20.txt' ) ) as f:
		Zlist = [];
		Alist = [];
		Nlist = [];
		amlist = [];
		amerrlist = [];
		amexplist = [];
		keylist = [];
		n = 0;
		lastZ = 0;
		for line in f:
			n = n+1;
			sZ   = line[0:6];
			sN   = line[7:14];
			sA   = line[15:22];
			Sym  = line[23:30];
			sam  = line[31:47];	# sam = standard atomic mass
			sam_err = line[48:56];
        
			# below is if there is a gap in one of the lines somewhere, I dont think it is needed for the new mass file?
			try:
				Z = int(sZ);
			except ValueError:
				Z = lastZ;
				pass;
        
			try:
				A = int(sA);
			except ValueError:
				A = 0;
				pass;
			try:
				N = int(sN);
			except ValueError:
				N = 0;
				pass;
            
			am,am_err,am_experimental,am_exp_err = decompose_atomic_mass( sam,sam_err );

			keylist.append(Z*1000+A);
			Zlist.append(Z);
			Nlist.append(N);
			Alist.append(A);
			amlist.append( am );
			amerrlist.append( am_err );
			amexplist.append( am_experimental );

			lastZ = Z;

        
		mass_table = np.empty( len(Zlist), dtype=[('key','u4'),('Z','u4'),('N','u4'),('A','u4'),('m','f8'),('m_err','f8'),('m_is_experimental','?') ] );
		mass_table['key'] = keylist;
		mass_table['Z'] = Zlist;
		mass_table['N'] = Nlist;
		mass_table['A'] = Alist;
		mass_table['m']     = amlist;
		mass_table['m_err'] = amerrlist;
		mass_table['m_is_experimental'] = amexplist;

		return mass_table;


mass_table = _initialize_module();

def get_m( A, Z ):
	#Returns atomic mass (value, error, experimental) for isotope A,Z.
	key = Z*1000+A;
	# index = np.searchsorted(mass_table['key'],key)
	index = np.argwhere(mass_table['key']==key).flatten()[0];
	rec = mass_table[index];
	if rec['key'] != key:
		raise Exception('Mass for A={:d},Z={:d} is not found'.format(A,Z) );
	return (rec['m'],rec['m_err'], rec['m_is_experimental'] );


#==============
#==============







def Where_things_done(angle_low, step_size, num_steps):
	# MAKE EDITS SO THAT NEUTRONS CAN BE ADDED AS A POSSIBLE BEAM OR TARGET THING

	# if Input_beam_Z.get() == str:	# just an idea to allow for symbol or Z inputting





	# the main inputs
	beam_A		=	int(Input_beam_A.get());
	beam_Z		=	int(Input_beam_Z.get());
	beam_energy	=	float(Input_beam_energy.get());		# MeV
	target_A	=	int(Input_target_A.get());
	target_Z	=	int(Input_target_Z.get());
	ejectile_A	=	int(Input_ejectile_A.get());
	ejectile_Z	=	int(Input_ejectile_Z.get());
	ejectile_ex_en	=	float(Input_ejectile_ex_en.get());	# MeV
	recoil_ex_en	=	float(Input_recoil_ex_en.get());	# MeV

	# working out the recoil stuff
	recoil_Z	=	int(beam_Z + target_Z - ejectile_Z);	# rZ = bZ+tZ-eZ
	recoil_A	=	int(beam_A + target_A - ejectile_A);	# rA = bA+tA-eA

	# names of the nuclei, uses IsotopeInfo.py to read masses.txt
	# beam_nucl	=	IsotopeInfo.get_sym(beam_Z);
	# target_nucl	=	IsotopeInfo.get_sym(target_Z);
	# ejectile_nucl	=	IsotopeInfo.get_sym(ejectile_Z);
	# recoil_nucl	=	IsotopeInfo.get_sym(recoil_Z);
	beam_nucl	=	get_sym(beam_Z);
	target_nucl	=	get_sym(target_Z);
	ejectile_nucl	=	get_sym(ejectile_Z);
	recoil_nucl	=	get_sym(recoil_Z);

	# Theta bits
	# angle_low	=	float(Input_angle_low.get());	# lab angle to start on
	# step_size	=	float(Input_step_size.get());
	# num_steps	=	int(Input_num_steps.get());
	angle_high	=	angle_low + ((num_steps+1)*step_size);
	lab_angle	=	np.arange(angle_low,angle_high,step_size);





	# Obtain atomic masses in amu and MeV
	beam_Mass_amu		= get_m(beam_A,beam_Z)[0];
	# beam_Mass_MeV		= IsotopeInfo.get_m(beam_A,beam_Z)[0]*physc.u_to_MeV;
	beam_Mass_MeV		= get_m(beam_A,beam_Z)[0]*u_to_MeV;
	# beam_DeltaM_MeV		= (beam_Mass_amu-beam_A)*physc.u_to_MeV;
	beam_DeltaM_MeV		= (beam_Mass_amu-beam_A)*u_to_MeV;

	target_Mass_amu		= get_m(target_A,target_Z)[0];
	# target_Mass_MeV		= IsotopeInfo.get_m(target_A,target_Z)[0]*physc.u_to_MeV;
	target_Mass_MeV		= get_m(target_A,target_Z)[0]*u_to_MeV;
	# target_DeltaM_MeV	= (target_Mass_amu-target_A)*physc.u_to_MeV;
	target_DeltaM_MeV	= (target_Mass_amu-target_A)*u_to_MeV;

	ejectile_Mass_amu	= get_m(ejectile_A,ejectile_Z)[0];
	# ejectile_Mass_MeV	= IsotopeInfo.get_m(ejectile_A,ejectile_Z)[0]*physc.u_to_MeV;
	ejectile_Mass_MeV	= get_m(ejectile_A,ejectile_Z)[0]*u_to_MeV;
	# ejectile_DeltaM_MeV	= (ejectile_Mass_amu-ejectile_A)*physc.u_to_MeV;
	ejectile_DeltaM_MeV	= (ejectile_Mass_amu-ejectile_A)*u_to_MeV;

	recoil_Mass_amu		= get_m(recoil_A,recoil_Z)[0];
	# recoil_Mass_MeV		= IsotopeInfo.get_m(recoil_A,recoil_Z)[0]*physc.u_to_MeV;
	recoil_Mass_MeV		= get_m(recoil_A,recoil_Z)[0]*u_to_MeV;
	# recoil_DeltaM_MeV	= (recoil_Mass_amu-recoil_A)*physc.u_to_MeV;
	recoil_DeltaM_MeV	= (recoil_Mass_amu-recoil_A)*u_to_MeV;

	# Q-value calculation
	GndSt_Qvalue	= beam_Mass_MeV + target_Mass_MeV - ejectile_Mass_MeV - recoil_Mass_MeV;	# Ground state Qvalue
	FinSt_Qvalue	= GndSt_Qvalue - ejectile_ex_en - recoil_ex_en;					# Final state Qvalue
	

	# y		= sum of masses of beam and target nuclei (MeV)
	# betac		= ?	something with beam energy and the masses of beam and target in MeV
	# ecmi		= ?	something with sum of masses (MeV) and 2(beam energy + target mass sum)
	# ecmf		= ?	sum of ecmi, final state Q value, ejectile mass (MeV), recoil mass (MeV) minus the sum of beam and target masses (MeV)
	# e3cm		= ?	uses ecmf, and ejectile and recoil masses (MeV)
	# beta3c	= ?	uses e3cm and mass of ejectile (MeV)
	# y_new		= ?	update to y?? uses e3cm, ejectile mass (MeV) and betac
	# cosagl	= cosine of the lab angle (need to translate angle to radians)
	# b		= term in bigger eqn?	negative betac multiplied by cosagl
	# a		= term in bigger eqn?	it is just y_new + b^2
	# c		= term in bigger eqn?	1 - y_new
	# d**2		= term in bigger eqn?	b^2 - (a*c)
	# b3L1		= ?	uses b, d**2, and a
	# y_loop	= ?	uses lots of above
	# y_adj		= y_loop = cos(COM angle)????
	# angl4		= recoil angle?
	# b3L2		=
	# y_loop2	=
	# y_adj2	=
	# angl4_2	=

	why	= beam_Mass_MeV + target_Mass_MeV;				# formerly "y"
	betac	= np.sqrt(beam_energy * (beam_energy + (2*beam_Mass_MeV)))/(why + beam_energy);
	ecmi	= np.sqrt( np.power(why,2) + (2* beam_energy * target_Mass_MeV) );
	ecmf	= ecmi + FinSt_Qvalue - why + ejectile_Mass_MeV + recoil_Mass_MeV;
	e3cm	= (np.power(ecmf,2)+((ejectile_Mass_MeV + recoil_Mass_MeV)*(ejectile_Mass_MeV - recoil_Mass_MeV)))/(2*ecmf);
	beta3c	= np.sqrt(1-(1/np.power(e3cm/ejectile_Mass_MeV,2)));
	why_new	= np.power(e3cm/ejectile_Mass_MeV,2)*(1-np.power(betac,2));	# formerly "why_new"

	cosagl	= np.cos(np.radians(lab_angle));
	bbb	= -betac*cosagl;		# formerly "b"
	aaa	= why_new + np.power(bbb,2);	# formerly "a"
	ccc	= 1 - why_new;			# formerly "c"
	ddd2	= np.power(bbb,2) - (aaa*ccc);	# formerly "d2" or "d**2"

	b3L1	= (-bbb + np.sqrt(ddd2))/aaa;
	# print(np.isnan(b3L1));		# I think it is somethink like "if b3L1 not a number then print that and set to -100" (returns array of indices where number = nan)
	b3L1[np.isnan(b3L1)]	= -100;	# setting to -100 if it is not a number ???

	lab_energy 	= ejectile_Mass_MeV*((1/np.sqrt(1-np.power(b3L1,2)))-1);
	lab_energy[b3L1<=0.0000001]	= np.nan;	# checking if b3L1 us less than 0.0000001, if so: not a number (similar to the b3L1 check)

	recoil_energy	= beam_energy + FinSt_Qvalue - lab_energy;	# formerly "rec_energy"

	why_loop	= (b3L1*cosagl-betac)/((1-(betac*b3L1*cosagl))*beta3c);	# formerly "y_loop"

	why_adj		= why_loop;
	why_adj[(np.abs(why_adj)>1) & (np.abs(why_adj)<=1.00001)]	= np.sign(why_adj[(np.abs(why_adj)>1) & (np.abs(why_adj)<=1.00001)]);
	why_adj[(np.abs(why_adj)>1.00001)]	= 0;

	cm_angle	= np.degrees(np.arccos(why_adj));	# Centre Of Mass Angle
	cm_angle[np.abs(why_adj+1)<0.000001]	= 180;
	cm_angle[np.abs(why_adj-1)<0.000001]	= 0;
	cm_angle[np.isnan(lab_energy)]		= np.nan;

	cm_lab_ratio	= beta3c * (1+(bbb/b3L1))/(1+(betac*beta3c*why_adj))/b3L1;
	cm_lab_ratio[np.isnan(lab_energy)]	= np.nan;

	# not sure what this is for? - must be important
	angl4	= np.degrees(np.arcsin(np.sqrt((lab_energy*(lab_energy+(2*ejectile_Mass_MeV)))/(recoil_energy*(recoil_energy+(2*recoil_Mass_MeV)))) * np.sin(np.radians(lab_angle)) ));	# IF YOU SEE THIS IN TERMINAL, ANSWER = NaN.

	recoil_angle	= angl4;	# formerly "rec_angle"
	check1	= lab_energy*(lab_energy+(2*ejectile_Mass_MeV)*np.power(cosagl,2))>(beam_energy*(beam_energy+(2*beam_Mass_MeV)));
	check2	= lab_angle<90;
	recoil_angle[(check1) & (check2)]	= 180-angl4[(check1) & (check2)];
	recoil_angle[np.isnan(lab_energy)]	= np.nan;

	# line 157 catkinpy.py
	output_array = np.array([lab_angle,cm_angle,lab_energy,cm_lab_ratio,recoil_angle,recoil_energy]);	# formerly "array"
	column_values = ['Lab angle','c.m angle','Lab energy','cm/lab ratio','Recoil angle','Recoil energy'];
	data_output = {'Lab angle':lab_angle,'c.m angle':cm_angle,'Lab energy':lab_energy,'cm/lab ratio':cm_lab_ratio,'Recoil angle':recoil_angle,'Recoil energy':recoil_energy};	# formerly "d"
	data_output_final = pd.DataFrame(data = data_output);	# formerly "df"





	# Returns something called a tuple	-	maybe need to set up as vector or big array when translating to other languages?
	return beam_A, beam_Z, beam_energy, beam_nucl, beam_Mass_MeV, target_A, target_Z, target_nucl, target_Mass_MeV, ejectile_A, ejectile_Z, ejectile_nucl, ejectile_Mass_MeV, recoil_A, recoil_Z, recoil_nucl, recoil_Mass_MeV, ejectile_ex_en, recoil_ex_en, GndSt_Qvalue, FinSt_Qvalue, data_output_final, lab_angle,cm_angle,lab_energy,cm_lab_ratio,recoil_angle,recoil_energy, why, betac, ecmi, ecmf, e3cm, beta3c, why_new;



def Director():

	# MAKE EDITS SO THAT NEUTRONS CAN BE ADDED AS A POSSIBLE BEAM OR TARGET THING

	# if Input_beam_Z.get() == str:	# just an idea to allow for symbol or Z inputting

	angle_low	=	float(Input_angle_low.get());	# lab angle to start on
	step_size	=	float(Input_step_size.get());
	num_steps	=	int(Input_num_steps.get());

	beam_A, beam_Z, beam_energy, beam_nucl, beam_Mass_MeV, target_A, target_Z, target_nucl, target_Mass_MeV, ejectile_A, ejectile_Z, ejectile_nucl, ejectile_Mass_MeV, recoil_A, recoil_Z, recoil_nucl, recoil_Mass_MeV, ejectile_ex_en, recoil_ex_en, GndSt_Qvalue, FinSt_Qvalue, data_output_final, lab_angle,cm_angle,lab_energy,cm_lab_ratio,recoil_angle,recoil_energy, why, betac, ecmi, ecmf, e3cm, beta3c, why_new = Where_things_done(angle_low, step_size, num_steps);


	## Below is the bit that writes the output:

	line_1 = f"Beam:\t A = {beam_A}, Z = {beam_Z}, {beam_nucl}, Mass = {beam_Mass_MeV} \n";
	line_2 = f"Target:\t A = {target_A}, Z = {target_Z}, {target_nucl}, Mass = {target_Mass_MeV} \n";
	line_3 = f"Ejectile:\t A = {ejectile_A}, Z = {ejectile_Z}, {ejectile_nucl}, Mass = {ejectile_Mass_MeV} \n";
	line_4 = f"Recoil:\t A = {recoil_A}, Z = {recoil_nucl}, Mass = {recoil_Mass_MeV} \n";
	line_5 = f"Beam E = {beam_energy} MeV\nEjectile Ex En = {ejectile_ex_en} MeV, Recoil Ex En = {recoil_ex_en} MeV \n";
	line_6 = f"Reaction Q-Value: Gnd State = {GndSt_Qvalue}, Final State = {FinSt_Qvalue}";

	Output_data_repeat["text"] = line_1 + line_2 + line_3 + line_4 + line_5 + line_6;

	# main_useful_1 = f"Lab angle (deg) \t COM angle (deg) \t Lab energy (MeV) \t COM/Lab ratio \t Recoil angle (deg) \t Recoil energy (MeV)\n";
	# main_useful_2 = f"{lab_angle} \t {cm_angle} \t {lab_energy} \t {cm_lab_ratio} \t {recoil_angle} \t {recoil_energy}";

	# Output_data["text"] = "TESTING SOMETHING\n this too\n wertyuiopdfghjklxcvbnm,kwhdiwgefohwaibf";
	# Output_data["text"] = main_useful_1 + main_useful_2;
	Output_data["text"] = data_output_final;

	# Output_data["text"] = f"BEAM A = {beam_A}, Beam Z = {beam_Z}, Element = {IsotopeInfo.get_sym(beam_Z)}";

	# print("test this");	# works

	# Output_data["text"] = multiple_outs();




	## Below prints the same output info to the terminal:
	print('#'*20);
	print('CATKINPY Calculation');
	print();
	print('Beam: ',beam_nucl);
	print('Z = {:d}, A = {:d}, Mass = {:f}'.format(beam_Z,beam_A,beam_Mass_MeV));
	print('Beam Energy = {:f}'.format(beam_energy));
	print();
	print('Target: ',target_nucl);
	print('Z = {:d}, A = {:d}, Mass = {:f}'.format(target_Z,target_A,target_Mass_MeV));
	print();
	print('Ejectile: ',ejectile_nucl);
	print('Z = {:d}, A = {:d}, Mass = {:f}'.format(ejectile_Z,ejectile_A,ejectile_Mass_MeV));
	print('Ejectile excitation = {:f}'.format(ejectile_ex_en));
	print();
	print('Recoil: ',recoil_nucl);
	print('Z = {:d}, A = {:d}, Mass = {:f}'.format(recoil_Z,recoil_A,recoil_Mass_MeV));
	print('Recoil excitation = {:f}'.format(recoil_ex_en));
	print('\n\n');
	print('Reaction Q-value');
	print('Ground state = {:f}'.format(GndSt_Qvalue));
	print('Final state = {:f}'.format(FinSt_Qvalue));

	print("\n"*3);
	print("Showing some mid-calculation steps ...");
	print(f"y = {why}");
	print(f"ecmi = {ecmi}");
	print(f"ecmf = {ecmf}");
	print(f"e3cm = {e3cm}");
	print(f"beta3c = {beta3c}");
	print(f"y_new = {why_new}");


	print("\n"*5);
	print(data_output_final);

	return # just in case things dont work, this cuts it off?





def Save_to_txt(angle_low, step_size, num_steps):
# def Save_to_txt():
	now = datetime.datetime.now();
	now_time_title = now.strftime("%Y_%m_%d__%H_%M_%S");
	now_time = now.strftime("%Y-%m-%d %H:%M:%S");
	time_date = f"Time and Date: {now_time}\n";

	Output_title = "Catkin with Python\nMade By Reuben, Based on CatkinPy by Jacob which was in turn based on Catkin By Wilton\n\n\n";

	output_txt = open(f"Output_Catkin_{now_time_title}.txt", "w+");



	beam_A, beam_Z, beam_energy, beam_nucl, beam_Mass_MeV, target_A, target_Z, target_nucl, target_Mass_MeV, ejectile_A, ejectile_Z, ejectile_nucl, ejectile_Mass_MeV, recoil_A, recoil_Z, recoil_nucl, recoil_Mass_MeV, ejectile_ex_en, recoil_ex_en, GndSt_Qvalue, FinSt_Qvalue, data_output_final, lab_angle,cm_angle,lab_energy,cm_lab_ratio,recoil_angle,recoil_energy, why, betac, ecmi, ecmf, e3cm, beta3c, why_new = Where_things_done(angle_low, step_size, num_steps);


	line_1 = f"Beam:\t A = {beam_A}, Z = {beam_Z}, {beam_nucl}, Mass = {beam_Mass_MeV} \n";
	line_2 = f"Target:\t A = {target_A}, Z = {target_Z}, {target_nucl}, Mass = {target_Mass_MeV} \n";
	line_3 = f"Ejectile:\t A = {ejectile_A}, Z = {ejectile_Z}, {ejectile_nucl}, Mass = {ejectile_Mass_MeV} \n";
	line_4 = f"Recoil:\t A = {recoil_A}, Z = {recoil_nucl}, Mass = {recoil_Mass_MeV} \n";
	line_5 = f"Beam E = {beam_energy} MeV\nEjectile Ex En = {ejectile_ex_en} MeV, Recoil Ex En = {recoil_ex_en} MeV \n";
	line_6 = f"Reaction Q-Value: Gnd State = {GndSt_Qvalue}, Final State = {FinSt_Qvalue}";
	line_7 = f"Showing some mid-calculation steps ...\ny = {why}\nbetac = {betac}\necmi = {ecmi}\necmf = {ecmf}\ne3cm = {e3cm}\nbeta3c = {beta3c}\ny_new = {why_new}\n"
	breaker = "\n"*3;
	information = breaker + line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + breaker + line_7 + breaker;


	saved_data = Output_title + time_date + information# + data_output_final;

	Output_data_repeat["text"] = "SAVED AS";	# changing the shown output on the screen
	Output_data["text"] = f"Output_Catkin_{now_time_title}.txt";

	output_txt.write(saved_data);
	output_txt.write(str(data_output_final));
	output_txt.close();
	return




#==============
#==============

# ADD IN A BIT TO SAVE THE RESULT AS A txt FILE
#	Make a button that says "save" (maybe also an entry box to specify what to save as?)
#	on clicking the button it does exactly the same as the regular button, but has extra loop to save things
#	that is it.

def Save_to_txt1():
	angle_low = float(Input_angle_low.get());
	step_size = float(Input_step_size.get());
	num_steps = int(Input_num_steps.get())
	return Save_to_txt(angle_low, step_size, num_steps);


def Save_to_txt2():
	angle_low = 0.0;
	step_size = 1.0;
	num_steps = 180;
	return Save_to_txt(angle_low, step_size, num_steps);

#==============
#==============


def Plotting_something():

	## setting up plot_window
	global plot_window;	# needs this or it can't find it in other functions
	plot_window = tk.Tk();	# declares window
	plot_window.title("CatKin Plot");
	# Window start sizing
	plot_window.rowconfigure(0, minsize = 10, weight = 1);
	plot_window.columnconfigure(0, minsize = 10, weight = 1);
	# plot_window.configure(bg='lightblue');


	## Getting everything defined and calculated before starting:
	angle_low	=	0.0	#float(Input_angle_low.get());	# lab angle to start on
	step_size	=	1.0	#float(Input_step_size.get());
	num_steps	=	180.0	#int(Input_num_steps.get());
	beam_A, beam_Z, beam_energy, beam_nucl, beam_Mass_MeV, target_A, target_Z, target_nucl, target_Mass_MeV, ejectile_A, ejectile_Z, ejectile_nucl, ejectile_Mass_MeV, recoil_A, recoil_Z, recoil_nucl, recoil_Mass_MeV, ejectile_ex_en, recoil_ex_en, GndSt_Qvalue, FinSt_Qvalue, data_output_final, lab_angle,cm_angle,lab_energy,cm_lab_ratio,recoil_angle,recoil_energy, why, betac, ecmi, ecmf, e3cm, beta3c, why_new = Where_things_done(angle_low, step_size, num_steps);

	## getting the selection of what to plot
	plot_selection = str(Plot_Menu_variable.get());


	## assigning what was chosen to what will be plotted.
	if plot_selection == "select option":	# Just a little reminder
		Output_data_repeat["text"] = "You need to select an output plot!"
		Output_data["text"] = "Fucking do it!!!!!";	# TEST LINE
		plot_something();
	if plot_selection == "Test":	# Just a little reminder
		test_plot();
	if plot_selection == "Ejectile - Energy vs Angle":
		ej_Lab_E_v_Lab_Ang(lab_energy, lab_angle);
	if plot_selection == "Recoil Lab Angle vs Ejectile Lab Angle":
		rec_Ang_v_Lab_Ang(recoil_angle, lab_angle);
	if plot_selection == "Recoil Energy vs Ejectile Lab Angle":
		rec_E_v_ej_Lab_Ang(recoil_energy, lab_angle);
	if plot_selection == "Recoil - Energy vs Angle":
		rec_E_v_rec_Ang(recoil_energy, recoil_angle);
	if plot_selection == "Lab angle vs COM Angle":
		Lab_Ang_v_COM_Ang(lab_angle, cm_angle, recoil_angle);
	if plot_selection == "K.E. vs COM Angle":
		KE_v_COM_Ang(lab_energy, cm_angle, recoil_energy);
	if plot_selection == "Lab Energy vs Lab Angle (Plectrum Plots)":
		Plectrum_Plots(lab_angle, lab_energy, recoil_angle, recoil_energy);


	# make the plot_window appear on the screen
	plot_window.mainloop();
	return;





def test_plot():
	fig1 = Figure(figsize = (5,5), dpi = 300); # makes a "figure"
	y = [i+1 for i in range(101)];	# defining y
	plot1 = fig1.add_subplot(111);	# placing plot on fig?
	plot1.plot(y);	#plot y
	canvas = FigureCanvasTkAgg(fig1,master=plot_window);	# assigning a canvas to tkinter
	canvas.draw();
	# canvas.get_tk_widget().pack();	#gluing the canvas with the stuff on onto the window
	toolbar = NavigationToolbar2Tk(canvas, plot_window);	#adding toolbar to the plot
	toolbar.update();
	canvas.get_tk_widget().pack();	#gluing on again after update
	return;

def plot_something():
	fig1 = Figure(figsize = (5,5), dpi = 300); # makes a "figure"
	y = np.sin([i for i in range(101)]);	# defining y
	plot1 = fig1.add_subplot(111);	# placing plot on fig?
	plot1.plot(y);	#plot y
	plot1.set_title("Just Select A Plot Option, Stupid.",fontsize=10);
	canvas = FigureCanvasTkAgg(fig1,master=plot_window);	# assigning a canvas to tkinter
	canvas.draw();
	# canvas.get_tk_widget().pack();	#gluing the canvas with the stuff on onto the window
	toolbar = NavigationToolbar2Tk(canvas, plot_window);	#adding toolbar to the plot
	toolbar.update();
	canvas.get_tk_widget().pack();	#gluing on again after update
	return;

# small plots:
#	ejectile lab energy vs lab angle	+ "second solution"
#	recoil angle vs recoil energy (says ejectile lab angle in CatKin)	+ "second solution"
#	recoil energy vs ejectile lab angle	+ "second solution"
#	recoil energy vs recoil angle		+ "second solution"
# Lab angle vs COM angle (ejectile & recoil)
# Lab K.E. vs COM angle (ejectile & recoil)
# "Plectrum Plots" - Lab angle vs Lab energy (ejectile & recoil)
# 

def ej_Lab_E_v_Lab_Ang(lab_energy, lab_angle):	# ejectile lab energy vs lab angle
	fig = Figure(figsize = (5,5), dpi = 300); # makes a "figure"
	y_values = lab_energy;	# defining y axis values
	x_values = lab_angle;
	plt = fig.add_subplot();
	plt.plot(x_values, y_values);
	plt.set_ylabel('ejectile lab energy [MeV]',fontsize=10);
	plt.set_xlabel('ejectile lab angle [deg]',fontsize=10, labelpad=0.5);	# labelpad needed to move x label up a bit
	plt.set_title("ejectile lab energy vs ejectile lab angle",fontsize=10);
	canvas = FigureCanvasTkAgg(fig,master=plot_window);
	canvas.draw();
	toolbar = NavigationToolbar2Tk(canvas, plot_window);
	toolbar.update();
	canvas.get_tk_widget().pack();
	return;

def rec_Ang_v_Lab_Ang(recoil_angle, lab_angle):	# recoil energy vs ejectile lab angle
	fig = Figure(figsize = (5,5), dpi = 300);
	y_values = recoil_angle;
	x_values = lab_angle;
	plt = fig.add_subplot();
	plt.plot(x_values, y_values);
	plt.set_ylabel('recoil angle [deg]',fontsize=10);
	plt.set_xlabel('ejectile lab angle [deg]',fontsize=10, labelpad=0.5);	# labelpad needed to move x label up a bit
	plt.set_title("recoil angle vs ejectile lab angle",fontsize=10);
	canvas = FigureCanvasTkAgg(fig,master=plot_window);
	canvas.draw();
	toolbar = NavigationToolbar2Tk(canvas, plot_window);
	toolbar.update();
	canvas.get_tk_widget().pack();
	return;

def rec_E_v_ej_Lab_Ang(recoil_energy, lab_angle):	#recoil energy vs ejectile lab angle
	fig = Figure(figsize = (5,5), dpi = 300);
	y_values = recoil_energy;
	x_values = lab_angle;
	plt = fig.add_subplot();
	plt.plot(x_values, y_values);
	plt.set_ylabel('recoil energy [MeV]',fontsize=10);
	plt.set_xlabel('ejectile lab angle [deg]',fontsize=10, labelpad=0.5);	# labelpad needed to move x label up a bit
	plt.set_title("recoil energy vs ejectile lab angle",fontsize=10);
	canvas = FigureCanvasTkAgg(fig,master=plot_window);
	canvas.draw();
	toolbar = NavigationToolbar2Tk(canvas, plot_window);
	toolbar.update();
	canvas.get_tk_widget().pack();
	return;

def rec_E_v_rec_Ang(recoil_energy, recoil_angle):	# recoil energy vs recoil angle
	fig = Figure(figsize = (5,5), dpi = 300);
	y_values = recoil_energy;
	x_values = recoil_angle;
	plt = fig.add_subplot();
	plt.plot(x_values, y_values);
	plt.set_ylabel('recoil energy [MeV]',fontsize=10);
	plt.set_xlabel('recoil angle [deg]',fontsize=10, labelpad=0.5);	# labelpad needed to move x label up a bit
	plt.set_title("recoil energy vs recoil angle",fontsize=10);
	canvas = FigureCanvasTkAgg(fig,master=plot_window);
	canvas.draw();
	toolbar = NavigationToolbar2Tk(canvas, plot_window);
	toolbar.update();
	canvas.get_tk_widget().pack();
	return;

def Lab_Ang_v_COM_Ang(lab_angle, cm_angle, recoil_angle):	# Lab angle vs COM angle (ejectile & recoil)
	fig = Figure(figsize = (5,5), dpi = 300);
	y_values_ej = lab_angle;	# ejectile Lab
	x_values = cm_angle;	# COM (remember- BOTH EJECTILE AND RECOIL USE THIS)
	y_values_rec = recoil_angle;	# recoil Lab
	plt = fig.add_subplot();
	plt.plot(x_values, y_values_ej, color="blue", label="Ejectile");
	plt.plot(x_values, y_values_rec, color="red", label="Recoil");
	plt.set_ylabel('Lab Angle [deg]',fontsize=10);
	plt.set_xlabel('COM Angle [deg]',fontsize=10, labelpad=0.5);	# labelpad needed to move x label up a bit
	plt.set_title("Lab Angle vs COM Angle",fontsize=10);
	plt.legend();
	canvas = FigureCanvasTkAgg(fig,master=plot_window);
	canvas.draw();
	toolbar = NavigationToolbar2Tk(canvas, plot_window);
	toolbar.update();
	canvas.get_tk_widget().pack();
	return;

def KE_v_COM_Ang(lab_energy, cm_angle, recoil_energy):	# Lab K.E. vs COM angle (ejectile & recoil)
	fig = Figure(figsize = (5,5), dpi = 300);
	y_values_ej = lab_energy;	# ejectile Lab energy
	x_values = cm_angle;	# COM (remember- BOTH EJECTILE AND RECOIL USE THIS)
	y_values_rec = recoil_energy;	# recoil Lab energy
	plt = fig.add_subplot();
	plt.plot(x_values, y_values_ej, color="blue", label="Ejectile");
	plt.plot(x_values, y_values_rec, color="red", label="Recoil");
	plt.set_ylabel('Lab K.E. [MeV]',fontsize=10);
	plt.set_xlabel('COM Angle [deg]',fontsize=10, labelpad=0.5);	# labelpad needed to move x label up a bit
	plt.set_title("Lab K.E. vs COM Angle",fontsize=10);
	plt.legend();
	canvas = FigureCanvasTkAgg(fig,master=plot_window);
	canvas.draw();
	toolbar = NavigationToolbar2Tk(canvas, plot_window);
	toolbar.update();
	canvas.get_tk_widget().pack();
	return;

def Plectrum_Plots(lab_angle, lab_energy, recoil_angle, recoil_energy):	# "Plectrum Plots" - Lab angle vs Lab energy (ejectile & recoil)
	fig = Figure(figsize = (5,5), dpi = 300);
	y_values_ej = lab_energy;	# ejectile Lab energy
	x_values_ej = lab_angle;	# ejectile Lab angle
	y_values_rec = recoil_energy;	# recoil Lab energy
	x_values_rec = recoil_angle;	# recoil Lab angle
	plt = fig.add_subplot();
	plt.plot(x_values_ej, y_values_ej, color="blue", label="Ejectile");
	plt.plot(x_values_rec, y_values_rec, color="red", label="Recoil");
	plt.set_ylabel('Lab Energy [MeV]',fontsize=10);
	plt.set_xlabel('Lab Angle [deg]',fontsize=10, labelpad=0.5);	# labelpad needed to move x label up a bit
	plt.set_title("Lab Energy vs Lab Angle",fontsize=10);
	plt.legend();
	canvas = FigureCanvasTkAgg(fig,master=plot_window);
	canvas.draw();
	toolbar = NavigationToolbar2Tk(canvas, plot_window);
	toolbar.update();
	canvas.get_tk_widget().pack();
	return;


#==============
#==============




#############################################################################################################
#############################################################################################################
#############################################################################################################

###############
# Windows & Buttons:
###############
window = tk.Tk();	# declares window
window.title("CatKin");

# Window start sizing
window.rowconfigure(0, minsize = 75, weight = 1);
window.columnconfigure(0, minsize = 50, weight = 1);
window.configure(bg='lightblue');
# window.configure(bg='white');
# window.configure(bg='lightpink');


#==============

TOP = tk.Frame(relief = tk.RAISED, master = window, bg="white");
INPUTS = tk.Frame(relief = tk.RAISED, master = window);		# Input values
BUTTON = tk.Frame(relief = tk.RAISED, master = window);		# Go button
OUTPUTS = tk.Frame(relief = tk.RAISED, master = window);	# Output Data
PLOTS = tk.Frame(relief = tk.RAISED, master = window, bg="lightblue");		# Where the Plots go
BOTTOM = tk.Frame(relief = tk.RAISED, master = window);		# bottom blank space, maybe good for credits

#==============

width1 = int(25);
width2 = int(5);
width3 = int(75);
total_width = int(150);

#==============


title_words = str("Gooey Cat");
words = str("Catkin kinematics calculator with a GUI [uses AME2020]\nCatKin [v2.03] (by Wilton) translated into python (by Jacob) then adapted and modified by me (Reuben).\n");
words2 = str("Use at your own peril, I assume this will work as needed ...");
words3 = str("If NaN comes up, it means not possible at that angle (ish).\nIf NaN for all: try swapping the input ejectile masses for the recoil masses.");
Top_Writing_words1 = tk.Label(text = title_words, master = TOP, bg = "white", width = int(total_width/2), height = 2, font = ("arial",15,"bold"));
Top_Writing_words2 = tk.Label(text = words+words2+words3, master = TOP, bg = "white", width = total_width, height = 5);
Top_Writing_words1.pack();
Top_Writing_words2.pack();


#==============

Input_beam_A_label	= tk.Label(master = INPUTS, text = 'A (beam)');	# yellow
Input_beam_A		= tk.Entry(master = INPUTS, width = width1);
Input_beam_Z_label	= tk.Label(master = INPUTS, text = 'Z (beam)');	# yellow
Input_beam_Z		= tk.Entry(master = INPUTS, width = width1);

Input_beam_energy_label	= tk.Label(master = INPUTS, text = 'Beam Energy [MeV]');	# green(light)
Input_beam_energy	= tk.Entry(master = INPUTS, width = width1);

Input_target_A_label	= tk.Label(master = INPUTS, text = 'A (target)');	# yellow
Input_target_A		= tk.Entry(master = INPUTS, width = width1);
Input_target_Z_label	= tk.Label(master = INPUTS, text = 'Z (target)');	# yellow
Input_target_Z		= tk.Entry(master = INPUTS, width = width1);

Input_ejectile_A_label	= tk.Label(master = INPUTS, text = 'A (ejectile)');	# yellow
Input_ejectile_A	= tk.Entry(master = INPUTS, width = width1);
Input_ejectile_Z_label	= tk.Label(master = INPUTS, text = 'Z (ejectile)');	#yellow
Input_ejectile_Z	= tk.Entry(master = INPUTS, width = width1);

Input_ejectile_ex_en_label	= tk.Label(master = INPUTS, text = 'Excitation Energy (ejectile) [MeV]');	# yellow
Input_ejectile_ex_en		= tk.Entry(master = INPUTS, width = width1);
Input_recoil_ex_en_label	= tk.Label(master = INPUTS, text = 'Excitation Energy (recoil) [MeV]');	# yellow
Input_recoil_ex_en		= tk.Entry(master = INPUTS, width = width1);

Input_angle_low_label	= tk.Label(master = INPUTS, text = 'Start Angle (lab) [deg]');	# yellow
Input_angle_low		= tk.Entry(master = INPUTS, width = width1);
Input_step_size_label	= tk.Label(master = INPUTS, text = 'Angle Step Size');	# yellow
Input_step_size		= tk.Entry(master = INPUTS, width = width1);
Input_num_steps_label	= tk.Label(master = INPUTS, text = 'Number of Steps');	# yellow
Input_num_steps		= tk.Entry(master = INPUTS, width = width1);



## a couple of bits to have a starting value and prevent errors
Input_beam_A.insert(0,"80");
Input_beam_Z.insert(0,"38");
Input_beam_energy.insert(0,"340");
# Input_target_A.insert(0,"196");	#	Pt target
# Input_target_Z.insert(0,"78");	#
Input_target_A.insert(0,"208");
Input_target_Z.insert(0,"82");
Input_ejectile_A.insert(0,"208");
Input_ejectile_Z.insert(0,"82");
Input_ejectile_ex_en.insert(0,"2.614");
Input_recoil_ex_en.insert(0,"0.386");
Input_angle_low.insert(0,"45");
Input_step_size.insert(0,"0.1");
Input_num_steps.insert(0,"10");







Input_beam_A_label.pack();
Input_beam_A.pack();
Input_beam_Z_label.pack();
Input_beam_Z.pack();
Input_beam_energy_label.pack();
Input_beam_energy.pack();
Input_target_A_label.pack();
Input_target_A.pack();
Input_target_Z_label.pack();
Input_target_Z.pack();
Input_ejectile_A_label.pack();
Input_ejectile_A.pack();
Input_ejectile_Z_label.pack();
Input_ejectile_Z.pack();
Input_ejectile_ex_en_label.pack();
Input_ejectile_ex_en.pack();
Input_recoil_ex_en_label.pack();
Input_recoil_ex_en.pack();
Input_angle_low_label.pack();
Input_angle_low.pack();
Input_step_size_label.pack();
Input_step_size.pack();
Input_num_steps_label.pack();
Input_num_steps.pack();





#==============

The_Button = tk.Button(text = "GO", font = ("default",20,"bold"), width = width2, height = 3, command = Director, master = BUTTON, bg = "green3");
The_Button.pack();

SAVE_Button = tk.Button(text = "Save to .txt", font = ("default",10), width = 10, height = 1, command = Save_to_txt1, master = BUTTON, bg = "blue");
SAVE_Button.pack();

SAVE_Button_FULL = tk.Button(text = "0-180 .txt", font = ("default",10), width = 10, height = 1, command = Save_to_txt2, master = BUTTON, bg = "red");
SAVE_Button_FULL.pack();

#==============

# repeating the key input information
Output_data_repeat = tk.Label(text = "--", fg = "black", bg = "white", width = width3, master = OUTPUTS);


# the useful output
Output_data_label = tk.Label(text = "Output", fg = "black", width = width3, master = OUTPUTS);			#, bg = "white"
Output_data = tk.Label(text = "--", fg = "black", bg = "white", width = width3, master = OUTPUTS);



Output_data_label.pack();
Output_data_repeat.pack();
Output_data.pack();


#==============

Plot_Menu_label = tk.Label(master = PLOTS, text = "Select Plot", width = width1);
# Plot_Menu_options = ["Test1", "Test2", "Test3"];
Plot_Menu_options = ["Test", "Ejectile - Energy vs Angle", "Recoil Lab Angle vs Ejectile Lab Angle", "Recoil Energy vs Ejectile Lab Angle", "Recoil - Energy vs Angle", "Lab angle vs COM Angle", "K.E. vs COM Angle", "Lab Energy vs Lab Angle (Plectrum Plots)"];
Plot_Menu_variable = tk.StringVar(PLOTS);
Plot_Menu_variable.set("select option");	# setting a starting label on the drop down menu
Plot_Menu = tk.OptionMenu(PLOTS, Plot_Menu_variable, *Plot_Menu_options);

Plot_Menu_label.pack();
Plot_Menu.pack();

Plot_Button = tk.Button(text = "Plot", font = ("default",10,"bold"), width = 10, height = 3, command = Plotting_something, master = PLOTS, bg = "green3");
Plot_Button.pack();


#==============
bottom_label = tk.Label(text = "Made by Reuben --> using Jacob's .py --> using Wilton's .xlsx", fg = "black", bg = "white", width = total_width, height = 1, master = BOTTOM);
bottom_label.pack();
#==============

TOP.grid(row = 0, column = 0, columnspan = 4);
INPUTS.grid(row = 1, column = 0);
BUTTON.grid(row = 1, column = 1);
OUTPUTS.grid(row = 1, column = 2);
PLOTS.grid(row = 1, column = 3);
BOTTOM.grid(row = 2, column = 0, columnspan = 4);

###############
# End:
###############
window.mainloop();	# opens window and runs the script on a loop until the window is closed.