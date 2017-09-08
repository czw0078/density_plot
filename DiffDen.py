#!/usr/bin/env python
import re
import os
import sys
import pickle
import math
import readline
import textwrap
import matplotlib
import matplotlib.mlab as mlab
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from function_read_file import *
from function_calculate_density import *
from optparse import OptionParser 
##########################################
# Analysis options                       #
##########################################
parser = OptionParser()	
parser.add_option('-x', dest = 'x_coords', help = 'x range in Amstrong', default = '-5.0 5.0 0.1')
parser.add_option('-y', dest = 'y_coords', help = 'y range in Amstrong', default = '-5.0 5.0 0.1')
parser.add_option('--p1', dest = 'p1', help = 'first point', default = False)
parser.add_option('--p2', dest = 'p2', help = 'second point', default = False)
parser.add_option('--p3', dest = 'p3', help = 'third point', default = False)
parser.add_option('-s', dest = 'show', help = 'label the values', default = False)
parser.add_option('-l', dest = 'levels', help = 'set total numbers of levels', default = False)
parser.add_option('-n', dest = 'nlines', help = 'numbers of lines.', default = "6")
parser.add_option('-k', dest = 'colour', help = 'grey or colour', default=False)
parser.add_option('-j', dest = 'judge', help = 'set criterian to pick points on plane', default=False)
parser.add_option('-d', dest = 'surface', help = 'genarate 3D surface instead of contour', default=False)
parser.add_option('-t', dest = 'tick', help = 'tick minor and major', default="0.5 1.0")
parser.add_option('-v', dest = 'q_range', help = 'range of the values', default=False)
parser.add_option('--heatmap', dest = 'heatmap', help = 'plot heatmap', default=False)
parser.add_option('--load', dest = 'load', help = 'load intermedia', default=False)
parser.add_option('--log', dest = 'log', help = 'logrithmatic', default=False)
params, args = parser.parse_args()
x_input=[ float(n) for n in params.x_coords.replace(',',' ').split() ]
y_input=[ float(n) for n in params.y_coords.replace(',',' ').split() ]
if params.p1==False or params.p2==False or params.p3==False:
	plane_setted=False
else:
	first_point=params.p1
	second_point=params.p2
	third_point=params.p3
	plane_setted=True

# -n analysis
try:
	tmp=params.nlines.split()
	if len(tmp)==1:
		n_line=int(params.nlines)
	else:
		n_plus,n_minus=[int(x) for x in tmp]
		n_line=n_plus+n_minus
except ValueError:
	n_line=6
# -l analysis
try:
	if params.levels!=False:
		v_limit=float(params.judge)
		nl=float(params.levels)
		delta2=v_limit/nl
		start2=-n_minus*delta2
		end2=n_plus*delta2
		i_levels=np.arange(start2,end2,delta2)	
except ValueError:
	params.levels=False
# -k analysis
is_colour=params.colour
# -j analysis
if params.judge != False:
	v_limit=float(params.judge)
else:
	v_limit=10.0
# -d analysis
is_surface=params.surface
# -t analysis
try:
	tick_minor,tick_major=params.tick.split()
	tick_minor=float(tick_minor)
	tick_major=float(tick_major)
except ValueError:
	tick_minor=0.5
	tick_major=1.0
# -v  analysis see params.q_range below
##########################################
# List objects                           #
##########################################
geom_info_list=[]
MO_info_list=[]
for each in args:
	tmp1,tmp2=auto_read(each)
	geom_info_list.append(tmp1)
	MO_info_list.append(tmp2)
n_obj=len(MO_info_list)
# test
# geom_info,MO_info=auto_read(args[0])
# print rho([0.0,1.0,0.0],MO_info)
# print geom_info
# exit()
# introduction
print "\n\n"+"-"*75+"\n\n"+" Differance Density Plot \n\n"+"-"*75
print textwrap.fill("The script now support both fchk format files from Gaussian \
and molden format from Molpro, Cfour. For post-hf method, \
make sure nature orbitals saved by using explicitly nature orbital keywords. \
Gaussian require another seperated job by keyword Guess to \
save nature orbitals into chk file before formatted. See workshop presentation.",75)
print "-"*75
print "%d objects found, and will plot the density as follow:\n\n\
%s \n " %(n_obj," - ".join(args))
print textwrap.fill("High-level calculation should be first therefore \
the subtraction is in correct order. To customized the density, \
change the order or use the old version of script.")

if not plane_setted:	
	# print coordinates
	l_format="{:>15}"*5
	print "\nMolecular coordinates as reference for setting coordinates of the plane:"
	print "-"*75
	print l_format.format("Name","Atomic","X","Y","Z")
	print "-"*75
	for each in geom_info_list[0]:
		print l_format.format(atomic2tag[each[1]],each[1],*each[0])
	print "-"*75+"\n"

##########################################
# Build mesh grid                        #
##########################################
if not plane_setted: # To do: put them into options
	# first point
	print textwrap.fill("First point on plane as origin point\
	of plane to be plot. Given multiple points, eg:(x1,y1,z1), (x2,y2,z2)...(xn,yn,zn), \
	the system will use the center point of them",75)
	first_point=raw_input("first >")
	# second point
	print ""
	print textwrap.fill("Input second point as x direction of plane. Multiple points input same as above",75)
	second_point=raw_input("second >")
	# third point
	print ""
	print textwrap.fill("Input third point on plane. Multiple points input same as above",75)
	third_point=raw_input("third >")

print "\nGenerating data..."
#set unit vectors for grid
p1,vx,vy=unit_vector( # bohr problem
	read_multi_xyz(first_point),
	read_multi_xyz(second_point),
	read_multi_xyz(third_point))
x_range=np.arange(*x_input)
y_range=np.arange(*y_input)
xx, yy=np.meshgrid(x_range, y_range)
v1,v2,v3=true_meshgrid(p1,vx,vy,xx,yy)
# see plane
# plt3d = plt.figure().gca(projection='3d')
# plt3d.plot_surface(v1, v2, v3)
# plt.show()
##########################################
# Generate Data                          #
##########################################
if params.load:
	# load from intermedia
	f_obj = open('intermedia.data')
	if f_obj:
		d1,d2,vv=pickle.load(f_obj)
	else:
		print "Can not load intermediate.data file"
		exit()
else:
	vv_list=[]
	for each_MO_info in MO_info_list:
		vv=rho([v1,v2,v3],each_MO_info)
		vv_list.append(vv)
	data_vv=vv_list[0]
	for each_vv in vv_list[1:]:
		data_vv=data_vv-each_vv
	# adjust
	d1=xx
	d2=yy
	vv=data_vv
	max_row=[]
	min_row=[]
	if params.q_range:
		for row in vv:
			max_row.append(max(row))
			min_row.append(min(row))
		print "Max:%s Min:%s \n"%(max(max_row),min(min_row))
		exit()
	# adjust value and cut too high number
	for i,row in enumerate(vv):
		for j,iterm in enumerate(row):
			if iterm>v_limit:
				vv[i][j]=v_limit
			elif iterm<-v_limit:
				vv[i][j]=-v_limit		
	data=[xx,yy,vv]
	print "Saving data to file intermedia.data"
	f_obj = open('intermedia.data', 'w')
	pickle.dump(data,f_obj)
	f_obj.close


##########################################
# Plot graph                             #
##########################################
fig=plt.figure(figsize=(11, 11))
plt.axis('scaled')
i_limit_x=max([abs(d1[0][0]),abs(d1[0][-1])])
i_limit_y=max([abs(d2[0][0]),abs(d2[-1][0])])
plt.axis([-i_limit_x,i_limit_x,-i_limit_y,i_limit_y])
# -l condition
if params.levels==False:
	last=n_line
else:
	last=i_levels
# --heatmap
if params.heatmap:
	plt.pcolor(d1, d2, vv, cmap='seismic',vmax=v_limit,vmin=-v_limit)#'RdBu'
	plt.colorbar()
	ax = fig.gca()
# 
elif is_surface:
	ax = fig.gca(projection='3d')
	surf=ax.plot_surface(d1, d2, vv, cmap=cm.coolwarm)
else:
	if is_colour: # -k condition
		CS = plt.contour(d1, d2, vv, last)#last here
	else:
		CS = plt.contour(d1, d2, vv, last, colors='k')
	ax = fig.gca()
	ax.set_xticks(np.arange(-i_limit_x, i_limit_x, tick_major))
	ax.set_xticks(np.arange(-i_limit_x, i_limit_x, tick_minor), minor = True)
	ax.set_yticks(np.arange(-i_limit_y, i_limit_y, tick_major))
	ax.set_yticks(np.arange(-i_limit_y, i_limit_y, tick_minor), minor = True)
	if params.show!=False:
		s=int(params.show)
		plt.clabel(CS,CS.levels[::s],inline=True,inline_spacing=0.01, fontsize=9,fmt='%.2e')
# Collect points needed to be display
print "Collect points needed to be display"
diplay_element_list=[]
diplay_data=[]
for each in geom_info_list[0]:
		tmp=each[0]-p1
		if on_plane(vx,vy,tmp,0.0001):
			a_xy=xy_new(vx,vy,tmp)
			if inside_limit(a_xy,i_limit_x,i_limit_y):
				tmp2=[a_xy,atomic2tag[each[1]]]
				diplay_data.append(a_xy)
				diplay_element_list.append(tmp2)

diplay_data=np.array(diplay_data)

# plt.scatter(diplay_data[:,0], diplay_data[:,1])
# diplay text

if not params.surface:
	for each in diplay_element_list:
		circle = plt.Circle((each[0][0], each[0][1]), 0.2, fc='w',lw=0)
		plt.gca().add_patch(circle)
		# 0.03 adjust of padding
		plt.text(each[0][0],each[0][1],each[1],
			horizontalalignment='center',verticalalignment='center')

# show graph
plt.show()
# example of usage:
#
# ./DiffDen.py A_B.molden A.molden B.molden 
# --p1 "0,0,0" --p2 "1,0,0" --p3 "1,0,1" -j 0.1 -d True -x "-6.0 6.0 0.05" -y "-6.0 6.0 0.05"
#
# ./DiffDen.py A_B.molden A.molden B.molden 
# --p1 "0,0,0" --p2 "1,0,0" --p3 "1,0,1" -j 0.1 -n 150 -t "0.1 0.5" -k True
#
#./DiffDen.py A_B.molden A.molden B.molden 
# --p1 "0.0,0.0,0.0" --p2 "1,0,0" --p3 "1,0,1" -j 0.01 -x "-6.0 6.0 0.03" -y "-6.0 6.0 0.03" --heatmap True
#
#./DiffDen.py A_B.molden A.molden B.molden --p1 "0.0,0.0,0.0" --p2 "1,0,0" --p3 "1,0,1" 
#-j 0.01 -x "-6.0 6.0 0.05" -y "-6.0 6.0 0.05" -n "9 4"  -l 20 --load True
#
#./DiffDen.py MOLDEN_NAT_CCSDT.molden MOLDEN_NAT_MP2.molden 
# -j 0.1 -x "-6.0 6.0 0.05" -y "-6.0 6.0 0.05" -d True
#
#./DiffDen.py MOLDEN_NAT_CCSDT.molden MOLDEN_NAT_MP2.molden -j 0.005 
#-x "-6.0 6.0 0.05" -y "-6.0 8.0 0.05" 
#--p1 "0.0 0.0 2.02022585638" --p2 "-1.37513850734 -1.74434977673 2.02022585638" 
#--p3 "0.381784167486 -5.28964990087 -1.70084505159" --heatmap True

