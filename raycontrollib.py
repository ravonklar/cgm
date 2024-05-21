import numpy as np
from scipy.stats import uniform
import yt
yt.enable_plugins()
import h5py
import tng_tools_portable_plus as tng
import trident
from unyt import unyt_array

halofilepath = '/home/ravonklar/setup/TNG100-1/cutouts/'

def galLoad(haloid):
	filename = f'{halofilepath}cutout_{haloid}.hdf5'
	ds = yt.load(filename)
	r_vir, primary_pos=ds.hsvals['Group_R_Crit200'].value, ds.hsvals['GroupPos'].value
	return ds, r_vir, primary_pos

def angleGenerator(points):
	#Generate fairly uniformly distributed points across a sphere.
	#points is the number of points on the equator of the sphere, and somewhat the meridians.
	ang = np.linspace(0, 2*np.pi, points)
	#The first angle that isn't 0 serves as a base unit
	angUnit = ang[1]
	#Determine the number of points that fit between 0 and pi
	levels = np.floor(np.pi/angUnit)
	#This step ensures an odd number, forcing one set of points to be on the equator
	levels = 2*np.ceil(levels/2)-1
	finalPoints = []
	#The condition is redundant, but ensures proper typing and values
	for level in range(np.floor(levels).astype(int)):
		#Determine the theta value for each level. Since the total is odd, angUnit/2
		#must be added in order to allow proper spacing
		theta = angUnit*level+angUnit/2
		#The number of points on each level is determined by the available space on it
		phi = np.linspace(0, 2*np.pi, np.floor(points*np.sin(theta)).astype(int))
		for x in phi[1:]:
						finalPoints.append((x, theta))
	#Returns the list containing the points in angular coordinates
	return finalPoints

def unitVectors(finalPoints):
	#Generate the cartesian unit vectors associated with the angular coordinates
	sep_unit_vectors=[]
	for entry in finalPoints:
		x, y, z=np.cos(entry[0])*np.sin(entry[1]), np.sin(entry[0])*np.sin(entry[1]), np.cos(entry[1])
		sep_unit_vectors.append(np.array([x, y, z]))
	sep_unit_vectors=np.array(sep_unit_vectors)
	return sep_unit_vectors

	#Deprecated
	#ips = r_vir*uniform.rvs(size=len(finalPoints))
	#secondAngles = np.pi*uniform.rvs(size=len(finalPoints))

def rayCalc(sep_unit_vector, r_vir):
	#Generate cartesian halo-centric coordinate ray from the given unit vector and virial radius
	#Generates an impact parameter and a second angle.
	#The coordinates so far are the direction to the closest approach of the ray to the halo center,
	#requiring an impact parameter to specify that distance. As there are an infinite number of rays
	#that go through that point at different angles, we specify one: the second angle.
	ip = r_vir*uniform.rvs()
	secAngle = np.pi*uniform.rvs()
	#Since previously we forced the number of points to be odd, we have no points exactly at 0 or pi theta.
	#As such, none of the unit vectors will be straight up in the z direction, so we are safe to take the
	#cross product between such a z vector and our unit vector to obtain a new vector along the xy plane,
	#which we normalize to remain as a unit vector.
	normalVec1 = np.cross(sep_unit_vector, [0,0,1])/np.linalg.norm(np.cross(sep_unit_vector, [0,0,1]))
	#We can take the cross product between our new unit vector and the original to receive a new unit vector
	#that forms a 2D plane tangent to the original ray when combined with the other unit vector.
	normalVec2 = np.cross(sep_unit_vector, normalVec1)/np.linalg.norm(np.cross(sep_unit_vector, normalVec1))
	#With this plane, we can generate a new vector on it using the second angle to give a ray direction
	#that's tangent to the impact point and also random.
	secondVec = np.cos(secAngle)*normalVec1 + np.sin(secAngle)*normalVec2
	#The ray goes out 8.2 virial radii in total, since that should far exceed the bounds of the halo
	startRay = ip*np.array(sep_unit_vector)-4.1*r_vir*np.array(secondVec)
	endRay = ip*np.array(sep_unit_vector)+4.1*r_vir*np.array(secondVec)
	return startRay, endRay
