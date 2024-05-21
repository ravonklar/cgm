import subprocess
import os
import sys
import numpy as np
import pickle

import pathlib
path = pathlib.Path().resolve()

halofilepath = '/home/ravonklar/setup/TNG100-1/cutouts/'
rayfilepath = '/scratch/ravonklar/'

halofiles = os.listdir(halofilepath)
haloids = {}
for filepath in halofiles:
	filepath = filepath.replace('cutout_', '')
	haloid = filepath.replace('.hdf5', '')
	startInd = 0

	rayfiles = os.listdir(rayfilepath)
	for filename in rayfiles:
		if filename.startswith(f'gal{haloid}'):
			filename = filename.replace('.hdf5', '')
			filename = filename.replace(f'gal{haloid}_', '')
			filename = int(filename)
			if filename > startInd:
				startInd = filename + 1

	haloids[haloid] = startInd


with open("inactiveHalos.pkl", "wb") as f:
	pickle.dump(haloids, f)
