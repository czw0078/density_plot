#!/usr/bin/env python
import re
import math
import os
import sys

def auto_read(obj):
	this_file=open(obj)
	if obj.split(".")[-1]=='fchk':
		geom_info,MO_info=read_gaussian_fchk(this_file)
	elif obj.split(".")[-1]=='molden':
		geom_info,MO_info=read_molden(this_file)
	else:
		print "unknow file format: ",obj, " exit"
		exit()
	return geom_info,MO_info
	#list_MO_coef,list_MO_occupy,list_AO = MO_info

def read_gaussian_fchk(this_file):
	##########################################
	# Expected angular momentun sequence     #
	##########################################
	i_ang_v_gaussian=[
		#S,P,D same
		[0,0,0],
		[1,0,0],[0,1,0],[0,0,1],
		[2,0,0],[0,2,0],[0,0,2],[1,1,0],[1,0,1],[0,1,1],
		#F different order 2
		[3,0,0],[0,3,0],[0,0,3],[1,2,0],[2,1,0],[2,0,1],
		[1,0,2],[0,1,2],[0,2,1],[1,1,1],
		#G different order 1
		[0,0,4],[0,1,3],[0,2,2],[0,3,1],[0,4,0],[1,0,3],
		[1,1,2],[1,2,1],[1,3,0],[2,0,2],[2,1,1],[2,2,0],
		[3,0,1],[3,1,0],[4,0,0],
		#H different order 1
		[0,0,5],[0,1,4],[0,2,3],[0,3,2],[0,4,1],[0,5,0],
		[1,0,4],[1,1,3],[1,2,2],[1,3,1],[1,4,0],[2,0,3],
		[2,1,2],[2,2,1],[2,3,0],[3,0,2],[3,0,2],[3,1,1],
		[3,2,0],[4,0,1],[4,1,0],[5,0,0]]
	# pure gaussian
	i_ang_v_pure_gaussian=[
		# 5d, NShell=-2
		[2,0],[2,1],[2,-1],[2,2],[2,-2],
		# 7f, NShell=-3
		[3,0],[3,1],[3,-1],[3,2],[3,-2],[3,3],[3,-3],
		# 10g, NShell=-4
		[4,0],[4,1],[4,-1],[4,2],[4,-2],[4,3],[4,-3],[4,4],[4,-4]]
	# convert from shell to list of anguler vector
	shell2v={
	0:[i_ang_v_gaussian[0]],
	1:i_ang_v_gaussian[1:4],
	2:i_ang_v_gaussian[4:10],
	3:i_ang_v_gaussian[10:20],
	4:i_ang_v_gaussian[20:36],
	-2:i_ang_v_pure_gaussian[0:5],
	-3:i_ang_v_pure_gaussian[5:12],
	-4:i_ang_v_pure_gaussian[12:21]}
	##########################################
	# gaussian format                        #
	##########################################
	atomic_mark=re.compile(r'^\s*Atomic\s*numbers\s*\S')
	coords_mark=re.compile(r'^\s*Current\s*cartesian\s*coordinates\s*\S')
	guess_mark=re.compile(r'^\s*Guess\s*=\s*\S')
	n_electrons_mark=re.compile(r'^\s*Number\s*of\s*electrons\s*\S')
	shell_type_mark=re.compile(r'^\s*Shell\s*types\s*\S')
	n_primitive_mark=re.compile(r'^\s*Number\s*of\s*primitives\sper\sshell\s*\S')
	any_letter_begin_mark=re.compile(r'^\s*[a-zA-Z]\s*\S')
	primitive_exponent_mark=re.compile(r'^\s*Primitive\s*exponents\s*\S')
	contraction_mark=re.compile(r'^\s*Contraction\s*coefficients\s*\S')
	contraction_sp_mark=re.compile(r'^\s*P\(S=P\)\s*Contraction\s*coefficients\s*\S')
	coordinate_mark=re.compile(r'^\s*Coordinates\s*of\s*each\s*shell\s*\S')
	MO_energy_mark=re.compile(r'^\s*Alpha\s*Orbital\s*Energies\s*\S')
	MO_coefficient_mark=re.compile(r'^\s*Alpha\s*MO\s*coefficients\s*\S')
	#from gaussian formatted checkpoint file we have shell list:
	#list_shell_fchk [ang_tag_no,xyz,primitive_parameters_list])
	#ang_tag_no 0=s, 1=p, -1=sp, 2=6d, -2=5d, 3=10f, -3=7f
	##########################################
	# Read line by line                      #
	##########################################
	#this_file=open(sys.argv[1])
	coords_list=[]
	atomic_list=[]
	shell_type_list=[]
	n_primitive_list=[]
	exponent_list=[]
	contraction_list=[]
	contraction_old_list=[]
	xyz_list=[]
	MO_energy_list=[]
	MO_coefficient_list=[]
	is_post_hf=False
	for i, this_line in enumerate(this_file):
		if guess_mark.match(this_line):# see weather it is post-hf
			is_post_hf=True  
		elif n_electrons_mark.match(this_line):# get n_electron
			n_electron=float(this_line.split()[-1])
			break
	# after electron found atomic number
	for i, this_line in enumerate(this_file):
		if atomic_mark.match(this_line):
			break
	# read atomic number
	for i, this_line in enumerate(this_file):
		if 	any_letter_begin_mark.match(this_line):
			break
		atomic_list=atomic_list+[int(x) for x in this_line.split()]
	# found coords
	for i, this_line in enumerate(this_file):
		if coords_mark.match(this_line):
			break
	# read coords
	for i, this_line in enumerate(this_file):
		if 	any_letter_begin_mark.match(this_line):
			break
		coords_list=coords_list+[float(x) for x in this_line.split()]
	# regroup
	coords_list=[coords_list[i:i+3] for i  in range(0, len(coords_list), 3)]
	# found shell_type
	for i, this_line in enumerate(this_file):
		if shell_type_mark.match(this_line):# goto shell type section
			break
	# Start from shell type section
	for i, this_line in enumerate(this_file):
		if n_primitive_mark.match(this_line): # meet n_primitive section
			break
		shell_type_list=shell_type_list+[int(x) for x in this_line.split()]
	# Start from n_primitive
	for i, this_line in enumerate(this_file):
		if any_letter_begin_mark.match(this_line):
			break
		n_primitive_list=n_primitive_list+[int(x) for x in this_line.split()]
	# Find primitive exponents
	for i, this_line in enumerate(this_file):
		if primitive_exponent_mark.match(this_line): # meet exponent section
			break
	# Start primitive exponents
	for i, this_line in enumerate(this_file):
		if contraction_mark.match(this_line):
			break
		exponent_list=exponent_list+[float(x) for x in this_line.split()]
	# then Contraction
	have_sp=False
	for i, this_line in enumerate(this_file):
		list_word=this_line.split()
		if contraction_sp_mark.match(this_line):
			contraction_old_list=contraction_list # save contraction
			contraction_list=[]
			have_sp=True
			continue
		if coordinate_mark.match(this_line): # meet coordinates section
			break
		contraction_list=contraction_list+[float(x) for x in this_line.split()]
	# then Coordinates
	for i, this_line in enumerate(this_file):
		if any_letter_begin_mark.match(this_line):
			break
		xyz_list=xyz_list+[float(x) for x in this_line.split()]
	# found energy first
	for i, this_line in enumerate(this_file):
		if MO_energy_mark.match(this_line):
			break
	# start energy
	for i, this_line in enumerate(this_file):
		if MO_coefficient_mark.match(this_line): # meet coefficeint section
			break
		MO_energy_list=MO_energy_list+[float(x) for x in this_line.split()]
	# then coefficeint
	for i, this_line in enumerate(this_file):
		if any_letter_begin_mark.match(this_line):
			break
		MO_coefficient_list=MO_coefficient_list+[float(x) for x in this_line.split()]
	##########################################
	# Rearrangement of data                  #
	##########################################
	# regroup each 3 coordintes into one set
	i=0;xyz_tmp_list=[];tmp3=[]
	for each in xyz_list:
		i=i+1
		tmp3.append(each)
		if i==3:
			xyz_tmp_list.append(tmp3)
			tmp3=[]
			i=0
	xyz_list=xyz_tmp_list
	# regroup the primitives
	new_contraction_list=[];new_exponent_list=[];new_contraction_sp_list=[]
	i_primitives=0
	for each in n_primitive_list:
		start=i_primitives
		i_primitives=i_primitives+each
		new_exponent_list.append(exponent_list[start:i_primitives])
		if have_sp:
			new_contraction_sp_list.append(contraction_list[start:i_primitives])
			new_contraction_list.append(contraction_old_list[start:i_primitives])
		else:
			new_contraction_list.append(contraction_list[start:i_primitives])
	contraction_list=new_contraction_list
	exponent_list=new_exponent_list
	contraction_sp_list=new_contraction_sp_list
	# make a AO_list
	AO_list=[]
	for i in range(len(shell_type_list)):
		if shell_type_list[i]==-1:
			for each_ang in shell2v[0]:
				AO_list.append(
					[each_ang,
					xyz_list[i],
					zip(exponent_list[i],contraction_list[i])])
			for each_ang in shell2v[1]:
				AO_list.append(
					[each_ang,
					xyz_list[i],
					zip(exponent_list[i],contraction_sp_list[i])])
		else:
			for each_ang in shell2v[shell_type_list[i]]:
				AO_list.append(
					[each_ang,
					xyz_list[i],
					zip(exponent_list[i],contraction_list[i])])
	# calculate the occupy number
	n_AO=len(AO_list)
	n_pair_electron=int(n_electron/2.0)
	if is_post_hf:
		n_occ_list=MO_energy_list
	else:
		n_occ_list=[2.0]*n_pair_electron+[0.0]*(n_AO-n_pair_electron)
	# regroup each n_AO coefficients
	i=0;MO_coefficient_list_new=[];tmp=[]
	for each in MO_coefficient_list:
		i=i+1
		tmp.append(each)
		if i==n_AO:
			MO_coefficient_list_new.append(tmp)
			tmp=[]
			i=0
	MO_coefficient_list=MO_coefficient_list_new
	this_file.close()
	geom_info=[[x,y] for x,y in zip(coords_list,atomic_list)]
	MO_info=(MO_coefficient_list,n_occ_list,AO_list)
	return geom_info,MO_info
	##########################################
	# unit test                              #
	##########################################
	# print len(shell_type_list)
	# print len(xyz_list)
	# print len(exponent_list)
	# print len(contraction_list)
	# print len(contraction_sp_list)
	# print len(MO_energy_list)
	# print len(MO_coefficient_list)
	# print len(AO_list)
	# print len(n_occ_list)
	# i=1
	# for each in AO_list:
	# 	print i,":",each
	# 	i=i+1

def read_molden(this_file):
	# cartesian
	i_ang_v_molden=[
		#S,P,D same
		[0,0,0],
		[1,0,0],[0,1,0],[0,0,1],
		[2,0,0],[0,2,0],[0,0,2],[1,1,0],[1,0,1],[0,1,1],
		#F different order 2
		[3,0,0],[0,3,0],[0,0,3],[1,2,0],[2,1,0],[2,0,1],
		[1,0,2],[0,1,2],[0,2,1],[1,1,1],
		#G different order 2
		[4,0,0],[0,4,0],[0,0,4],[3,1,0],[3,0,1],[1,3,0],
		[0,3,1],[1,0,3],[0,1,3],[2,2,0],[2,0,2],[0,2,2],
		[2,1,1],[1,2,1],[1,1,2]]
	# pure same as gaussian
	i_ang_v_pure_molden=[
		# 5d, NShell=-2
		[2,0],[2,1],[2,-1],[2,2],[2,-2],
		# 7f, NShell=-3
		[3,0],[3,1],[3,-1],[3,2],[3,-2],[3,3],[3,-3],
		# 10g, NShell=-4
		[4,0],[4,1],[4,-1],[4,2],[4,-2],[4,3],[4,-3],[4,4],[4,-4]]
	# from label to vecter
	c_shell2v={
	's':[i_ang_v_molden[0]],
	'p':i_ang_v_molden[1:4],
	'd':i_ang_v_molden[4:10],
	'f':i_ang_v_molden[10:20],
	'g':i_ang_v_molden[20:36]}
	##########################################
	# Molden format                          #
	##########################################
	# strs="0.1240000000D+01  0.0000000000D+00"
	# [float(x.replace("D","e")) for x in strs.split()]
	atom_mark=re.compile(r'^\s*\[\s*[Aa][Tt][Oo][Mm][Ss]\s*\]\s*')
	gto_mark=re.compile(r'^\s*\[\s*GTO\s*\]\s*$')
	MO_mark=re.compile(r'^\s*\[\s*MO\s*\]\s*$')
	any_mark=re.compile(r'^\s*\[\s*[A-Za-z]*\s*\]\s*$')
	any_mark2=re.compile(r'^\s*\[\s*')# because cfour will have [Molden Format] twice
	cfour_mark=re.compile(r'^\s*\[\s*[A-Za-z]*\s*\]\s*$')
	occup_line=re.compile(r'^\s*[Oo][Cc][Cc][Uu][Pp]=\s*\d\.\d+\s*')
	atom_line=re.compile(r'^\s*\d+\s*0\s*$')
	shell_line=re.compile(r'^\s*[A-Za-z]+\s*\d+\s*1.00\s*$')
	coefficient_line=re.compile(r'^\s*\d+\s*-?\d\.\d+\s*$')# here the 10 digits for molpro and 11 digits for cfour
	#
	d5_start_mark1=re.compile(r'^\s*\[\s*5D\s*\]\s*$')
	d5_start_mark2=re.compile(r'^\s*\[\s*5D\s*7F\s*\]\s*$')
	d5_only_mark=re.compile(r'^\s*\[\s*5D\s*10F\s*\]\s*$')
	#
	f7_start_mark1=re.compile(r'^\s*\[\s*7F\s*\]\s*$')
	f7_start_mark2=re.compile(r'^\s*\[\s*7F\s*9G\s*\]\s*$')
	f7_only_mark=re.compile(r'^\s*\[\s*7F\s*15G\s*\]\s*$')
	#
	g9_start_mark=re.compile(r'^\s*\[\s*9G\s*\]\s*$')
	##########################################
	# File IO                                #
	##########################################
	# this_file=open(sys.argv[1])
	# goto [atoms] section
	for i, this_line in enumerate(this_file):
		if atom_mark.match(this_line):
			break

	# read coordinates of atoms !!!!!!
	list_atom=[]
	angstrom2bohr=1.88973
	for i, this_line in enumerate(this_file):
		if any_mark2.match(this_line):
			break
		list_word=this_line.split()
		if len(list_word)!=0:
			xyz=[float(x)*angstrom2bohr for x in list_word[3:6]]
			atomic_number=int(list_word[2])
			atomic_tag=list_word[0]
			list_atom.append([xyz,atomic_number,atomic_tag])
	
	# found [GTO] section
	this_file.seek(0)
	for i, this_line in enumerate(this_file):
		if d5_start_mark1.match(this_line) or d5_start_mark2.match(this_line):
			c_shell2v['d']=i_ang_v_pure_molden[0:5]
			c_shell2v['f']=i_ang_v_pure_molden[5:12]
			c_shell2v['g']=i_ang_v_pure_molden[12:21]
		elif f7_start_mark1.match(this_line) or f7_start_mark1.match(this_line):
			c_shell2v['f']=i_ang_v_pure_molden[5:12]
			c_shell2v['g']=i_ang_v_pure_molden[12:21]
		elif g9_start_mark.match(this_line):
			c_shell2v['g']=i_ang_v_pure_molden[12:21]
		elif d5_only_mark.match(this_line):
			c_shell2v['d']=i_ang_v_pure_molden[0:5]
		elif f7_only_mark.match(this_line):
			c_shell2v['f']=i_ang_v_pure_molden[5:12]
		elif gto_mark.match(this_line):
			break
	# flattern nested loop by store in varibles
	a_shell=[]
	s_shell=[]# for s part of sp shell
	p_shell=[]# for p part of sp shell
	# group shell together, assigned angular momentum and center
	# coordinates to be transformed to basis functions.
	list_shell=[]
	is_SP_schell=False
	for i, this_line in enumerate(this_file):
		list_word=this_line.split()
		if any_mark.match(this_line):		
			break
		elif atom_line.match(this_line):
			i_atom=int(list_word[0])
		elif shell_line.match(this_line):
			c_shell=list_word[0]
			a_shell=[c_shell,i_atom]
			i_primitives=int(list_word[1])
			if c_shell=="sp" or c_shell=="SP":
				is_SP_schell=True
				s_shell=['s',i_atom]
				p_shell=['p',i_atom]
		elif len(list_word)!=0:
			i_primitives=i_primitives-1
			if is_SP_schell:
				exponent,contraction_s,contraction_p=\
				[float(x.replace("D","e")) for x in list_word]
				if contraction_s!=0:
					s_shell.append([exponent,contraction_s])
				if contraction_p!=0:
					p_shell.append([exponent,contraction_p])
				if i_primitives==0:
					list_shell.append(s_shell)
					list_shell.append(p_shell)
					s_shell=[]
					p_shell=[]
					is_SP_schell=False
			else:
				(exponent,contraction)=\
				[float(x.replace("D","e")) for x in list_word]
				if contraction!=0:
					a_shell.append([exponent,contraction])
				if i_primitives==0:
					list_shell.append(a_shell)
					a_shell=[]
		else:
			continue
	# change list of shell into list of AO
	list_AO=[]
	for each_shell in list_shell:
		ang_v=c_shell2v[each_shell[0]]
		for each_ang in ang_v:
			list_AO.append([each_ang,
				list_atom[each_shell[1]-1][0],each_shell[2:]]) #??
	# found [MO] section
	this_file.seek(0)
	for i, this_line in enumerate(this_file):
		if MO_mark.match(this_line):
			break
	# read MO coefficients
	list_MO_coef=[]
	list_MO_occupy=[]
	each_MO_coef=[]
	for i, this_line in enumerate(this_file):
		list_word=this_line.split()
		if occup_line.match(this_line):
			if len(each_MO_coef)!=0:
				list_MO_coef.append(each_MO_coef)
				each_MO_coef=[]
			list_MO_occupy.append(float(list_word[1]))
		elif coefficient_line.match(this_line):
			each_MO_coef.append(float(list_word[1]))
		else:
			continue
	list_MO_coef.append(each_MO_coef)
	this_file.close()
	MO_info=(list_MO_coef,list_MO_occupy,list_AO)
	return list_atom,MO_info

atomic2tag=["X",
	"H", 	"He",	"Li",	"Be",	"B",	"C",	"N",	"O",	"F",	"Ne",
	"Na",	"Mg",	"Al",	"Si",	"P",	"S",	"Cl",	"Ar",	"K",	"Ca",
	"Sc",	"Ti",	"V",	"Cr",	"Mn",	"Fe",	"Co",	"Ni",	"Cu",	"Zn",
	"Ga",	"Ge",	"As",	"Se",	"Br",	"Kr",	"Rb",	"Sr",	"Y",	"Zr",
	"Nb",	"Mo",	"Tc",	"Ru",	"Rh",	"Pd",	"Ag",	"Cd",	"In",	"Sn",
	"Sb",	"Te",	"I",	"Xe",	"Cs",	"Ba",	"La",	"Ce",	"Pr",	"Nd",
	"Pm",	"Sm",	"Eu",	"Gd",	"Tb",	"Dy",	"Ho",	"Er",	"Tm",	"Yb",
	"Lu",	"Hf",	"Ta",	"W",	"Re",	"Os",	"Ir",	"Pt",	"Au",	"Hg",
	"Tl",	"Pb",	"Bi",	"Po",	"At",	"Rn",	"Fr",	"Ra",	"Ac",	"Th",
	"Pa",	"U"	,	"Np",	"Pu",	"Am",	"Cm",	"Bk",	"Cf",	"Es",	"Fm",
	"Md",	"No",	"Lr",	"Rf",	"Db",	"Sg",	"Bh",	"Hs",	"Mt",	"Ds",
	"Rg",	"Cn",	"Uut",	"Fl",	"Uup",	"Lv",	"Uus",	"Uuo"]

