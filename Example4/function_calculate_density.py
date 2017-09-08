#!/usr/bin/env python
import re
import math
import numpy as np

# xy within the limits
def inside_limit(xy,x_limit,y_limit):
	x,y=xy
	if abs(x)<x_limit and abs(y)<y_limit:
		return True
	else:
		return False

# decomposed
def xy_new(vx,vy,v):
	new_x=np.dot(vx,v)
	new_y=np.dot(vy,v)
	return np.array([new_x,new_y])

# mixed product
def on_plane(vx,vy,v,delta):
	tmp=np.cross(vx,vy)
	tmp=np.dot(tmp,v)
	if abs(tmp)<delta:
		return True
	else:
		return False

# unit vectors for grid
def unit_vector(p1,p2,p3):
	vo=p1
	vx=p2-vo
	vx=vx/np.linalg.norm(vx)
	vy=p3-vo
	vy=vy-np.dot(vx,vy)*vx
	vy=vy/np.linalg.norm(vy)
	return (p1,vx,vy)
# get real 3d meshgrid from unit vectors and origin
def true_meshgrid(p1,vx,vy,xx,yy):
	v1=xx*vx[0]+yy*vy[0]+p1[0]
	v2=xx*vx[1]+yy*vy[1]+p1[1]
	v3=xx*vx[2]+yy*vy[2]+p1[2]
	return (v1,v2,v3)
# read coordinates of plane to be plot
def read_multi_xyz(prompt):
	tmp=[float(x) for x in re.findall(r'-?\d+\.?\d*',prompt)]
	tmp=[np.array(tmp[i:i+3]) for i  in range(0, len(tmp), 3)]
	return sum(tmp)/float(len(tmp))
# from normalized to normalized on need alpha 
# see multiwfn fileIO.f90 and Git HORTON project
def pure_polys(ang, xyz):
	x,y,z=xyz
	poly_dict={
		# d=2
		(2,0):-0.5*x*x-0.5*y*y+z*z,
		(2,1):x*z,
		(2,-1):y*z,
		(2,2):3.0**0.5/2.0*(x*x-y*y),
		(2,-2):x*y,
		# f=3
		(3,0):-3.0/(2.0*5.0**0.5)*(x*x*z+y*y*z)+z*z*z,
		(3,1):-(3.0/8.0)**0.5*x*x*x-(3.0/40.0)**0.5*x*y*y+(6.0/5.0)**0.5*x*z*z,
		(3,-1):-(3.0/40.0)**0.5*x*x*y-(3.0/8.0)**0.5*y*y*y+(6.0/5.0)**0.5*y*z*z,
		(3,2):3**0.5/2.0*(x*x*z-y*y*z),
		(3,-2):x*y*z,
		(3,3):(5.0/8.0)**0.5*x*x*x-3.0/8.0**0.5*x*y*y,
		(3,-3):3.0/8.0**0.5*x*x*y-(5.0/8.0)**0.5*y*y*y,
		# g=4
		(4,0):z*z*z*z+3.0/8.0*(x*x*x*x+y*y*y*y)-3.0*(3.0/35.0)**0.5*(x*x*z*z+y*y*z*z-1.0/4.0*x*x*y*y),
		(4,1):2.0*(5.0/14.0)**0.5*x*z*z*z-3.0/2.0*(5.0/14.0)**0.5*x*x*x*z-3.0/2.0/14.0**0.5*x*y*y*z,
		(4,-1):2.0*(5.0/14.0)**0.5*y*z*z*z-3.0/2.0*(5.0/14.0)**0.5*y*y*y*z-3.0/2.0/14.0**0.5*x*x*y*z,
		(4,2):3.0*(3.0/28.0)**0.5*(x*x*z*z-y*y*z*z)-5.0**0.5/4.0*(x*x*x*x-y*y*y*y),
		(4,-2):3.0/7.0**0.5*x*y*z*z-(5.0/28.0)**0.5*(x*x*x*y+x*y*y*y),
		(4,3):(5.0/8.0)**0.5*x*x*x*z-3.0/8.0**0.5*x*y*y*z,
		(4,-3):-(5.0/8.0)**0.5*y*y*y*z+3.0/8.0**0.5*x*x*y*z,
		(4,4):35.0**0.5/8.0*(x*x*x*x+y*y*y*y)-3.0/4.0*3.0**0.5*x*x*y*y,
		(4,-4):5.0**0.5/2.0*(x*x*x*y-x*y*y*y)
		}
	return poly_dict[tuple(ang)]
##########################################
# Form MO by AOs                         #
##########################################
# remember As applied to computers, 
# IEEE 754 recommends several functions for computing a power.
# It defines 0**0=1
# normalizatin of each primitive
def normal_prim(lx,ly,lz,alpha):
	tmp=(2*alpha/math.pi)**(3.0/4.0)*(8*alpha)**((lx+ly+lz)/2.0)
	tmp=tmp*(
		math.factorial(lx)*1.0/math.factorial(2*lx)*
		math.factorial(ly)*1.0/math.factorial(2*ly)*
		math.factorial(lz)*1.0/math.factorial(2*lz))**0.5
	return tmp
# must consider both pure and cartesian form
def calc_bf(xyz,each_AO):
		(x,y,z)=xyz
		ang,(Ox,Oy,Oz),list_primitives=each_AO
		(rx,ry,rz)=[x-Ox,y-Oy,z-Oz]
		rr=rx**2+ry**2+rz**2
		exponent_term=0
		if len(ang)==2:# spherical
			poly_term=pure_polys(ang,xyz)
			for each_primitive in list_primitives:
				(alpha,contraction)=each_primitive
				exponent_term=normal_prim(sum(ang),0,0,alpha)*contraction\
				*math.e**(-alpha*rr)+exponent_term
		else: # cartesian
			(lx,ly,lz)=ang
			poly_term=rx**lx*ry**ly*rz**lz
			for each_primitive in list_primitives:
				(alpha,contraction)=each_primitive
				exponent_term=normal_prim(lx,ly,lz,alpha)*contraction\
				*math.e**(-alpha*rr)+exponent_term
		val=exponent_term*poly_term
		# print each_AO # debug
		return val
# test bf pure and cartesian
# print calc_bf([0.0,1.0,0.0],AO_list[126])
# print calc_bf([0.0,1.0,0.0],AO_list[49])
# print calc_bf([0.0,1.0,0.0],AO_list[4])
def construct_MO(xyz,list_MO_coef,list_MO_occupy,list_AO):
	wfn=[]
	for each_col,occupy in zip(list_MO_coef,list_MO_occupy):
		tmp=0
		if occupy==0:
			wfn.append(0.0)
			continue
		for each_coef,each_AO in zip(each_col,list_AO):
			tmp=tmp+each_coef*calc_bf(xyz,each_AO)
		wfn.append(tmp)
	return wfn
# test construct_MO
# a_wfn=construct_MO([0.0,1.0,0.0],MO_coefficient_list,AO_list)
# i=0
# for each in a_wfn:
# 	i=i+1
# 	print i,":",each
# exit()
def rho(xyz,MO_info):
	val=0
	list_MO_coef,list_MO_occupy,list_AO=MO_info
	wfn=construct_MO(xyz,list_MO_coef,list_MO_occupy,list_AO)
	for each_wfn,each_occupy in zip(wfn,list_MO_occupy):
		val=val+each_wfn**2.0*each_occupy
	return val
#test the density 2.57567729086e-05 
#(unit of [0.0,1.0,0.0] is Bohr)
# print rho([0.0,1.0,0.0],MO_coefficient_list,n_occ_list,AO_list)
# exit()