import subprocess
import time

import pathlib
path = pathlib.Path().resolve()

print(path)

import raycontrollib as rclib

import sys
import pickle
import os

points = 127

uvPoints = rclib.angleGenerator(points)

uvList = rclib.unitVectors(uvPoints)

totraynum = len(uvList)

remain = totraynum
procNum = 5

activeProc = []
#haloids = pickle.load(open('done_halos84', 'rb'))
#print(haloids)

halofiles = os.listdir('./TNG100-1/cutouts/84/')
haloids = []
for filepath in halofiles:
	filepath = filepath.replace('halo_', '')
	filepath = filepath.replace('.hdf5', '')
	haloids.append(filepath)

timeList = []

for haloid in haloids:
	print(f'BEGAN HALO {haloid}')
	boolean = True
	remain = totraynum
	i = 0
	while boolean:
		targetProc = min(remain, procNum)
		while len(activeProc) != targetProc:
			seed = time.perf_counter()
			result = subprocess.Popen(f'cd {path} && python rccreate.py {seed} {int(haloid)} {uvList[i][0]} {uvList[i][1]} {uvList[i][2]} {i}', shell=True, universal_newlines=True, stdout=subprocess.DEVNULL)
			activeProc.append(result)
			print(f'STARTED RAY WITH DATA {seed} {int(haloid)} {uvList[i][0]} {uvList[i][1]} {uvList[i][2]} {i}')
			i += 1
			targetProc = min(remain, procNum)

			timeList.append(seed)

		for proc in activeProc:
			#print(f'READING OUTPUT FROM {activeProc.index(proc)}')
			#print(proc.communicate())
			test = proc.poll()
			if test != None:
				indP = activeProc.index(proc)
				endTime = time.perf_counter()
				print(f'RAY COMPLETED IN {endTime-timeList[indP]}')
				timeList.pop(indP)
				activeProc.remove(proc)
				remain-=1
		if remain == 0:
			boolean = False
			print('Mission success')
	print(f'FINISHED HALO {haloid}')
	#SINGLE USE, REMEMBER TO DISABLE
	break