import subprocess
import time

import pathlib
path = pathlib.Path().resolve()

print(path)

rayfilepath = '/scratch/ravonklar/'

import raycontrollib as rclib

import sys
import pickle
import os

import h5py 
import requests
import astropy
import numpy as np
import tng_tools_portable_plus as tng
import yt
import trident
yt.enable_plugins()
from unyt import unyt_array


#import shutil

haloid = int(sys.argv[1])

startInd = int(sys.argv[2])

increment = int(sys.argv[3])

points = 127
uvPoints = rclib.angleGenerator(points)
uvList = rclib.unitVectors(uvPoints)

totraynum = len(uvList)
if startInd+increment > totraynum:
	increment = totraynum-startInd
uvPoints = uvPoints[startInd:startInd+increment]
uvList = uvList[startInd:startInd+increment]

totraynum = len(uvList)

remain = totraynum
procNum = 1

'''tempFiles = []
for num in range(procNum):
	shutil.copyfile(f'C:/Users/ryker/Desktop/CGM/suffer/halo_{haloid}.hdf5', f'C:/Users/ryker/Desktop/CGM/suffer/halo_{haloid}_{num}.hdf5')
	tempFiles.append(f'C:/Users/ryker/Desktop/CGM/suffer/halo_{haloid}_{num}.hdf5')'''

activeProc = []
#haloids = pickle.load(open('done_halos84', 'rb'))
#print(haloids)

timeList = []

print(f'BEGAN HALO {haloid} INDEX {startInd}')

#more "physical" variables
omega_L=0.6911
omega_m=0.3089
omega_b=0.0486
h=0.6774
kpctokm=3.086e16
kpctocm=kpctokm*1e5
#uni_age=cosmo.age(0)

ds, r_vir, primary_pos = rclib.galLoad(haloid)

for i in range(totraynum):

	sep_u_v = [uvList[i][0], uvList[i][1], uvList[i][2]]
	rayid = i+startInd

	seed = time.perf_counter()
	np.random.seed(int(seed))
	startRay, endRay = rclib.rayCalc(sep_u_v, r_vir)

	start=unyt_array(startRay, 'code_length', registry=ds.unit_registry)
	end=unyt_array(endRay, 'code_length', registry=ds.unit_registry)
	ray_filename=f'gal{haloid}_{rayid}'
		#POSSIBLY IMPORTANT: It looks like passing our start and end points into trident.make_simple_ray actually *changes* what our start and end points are defined as. For this reason, I re-define them here.
	start=unyt_array(startRay, 'code_length', registry=ds.unit_registry)
	end=unyt_array(endRay, 'code_length', registry=ds.unit_registry)

	tng.ray_func(ds, start_position=start, end_position=end, sn_ratio=18, complete_filename=f'{rayfilepath}{ray_filename}.hdf5')
	print(f'RAY COMPLETED IN {time.perf_counter()-seed}')

print(f'FINISHED HALO {haloid} END {startInd+increment-1}')


'''for file in tempFiles:
	os.remove(file)'''