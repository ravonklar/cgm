import h5py 
import requests
import pickle
import astropy
import pickle
import numpy as np
import tng_tools_portable_plus as tng
import yt
import trident
yt.enable_plugins()
from unyt import unyt_array

import raycontrollib as rclib
import sys


#more "physical" variables
omega_L=0.6911
omega_m=0.3089
omega_b=0.0486
h=0.6774
kpctokm=3.086e16
kpctocm=kpctokm*1e5
#uni_age=cosmo.age(0)

seed = sys.argv[1]
haloid = sys.argv[2]
sep_u_v = [float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])]
rayid = sys.argv[6]

seed = seed.replace('\n', '')
seed = seed.replace('.', '')
np.random.seed(int(seed[:8]))

ds, r_vir, primary_pos = rclib.galLoad(haloid)

startRay, endRay = rclib.rayCalc(sep_u_v, r_vir)

start=unyt_array(startRay, 'code_length', registry=ds.unit_registry)
end=unyt_array(endRay, 'code_length', registry=ds.unit_registry)
ray_filename=f'gal{haloid}_{rayid}'
	#POSSIBLY IMPORTANT: It looks like passing our start and end points into trident.make_simple_ray actually *changes* what our start and end points are defined as. For this reason, I re-define them here.
start=unyt_array(startRay, 'code_length', registry=ds.unit_registry)
end=unyt_array(endRay, 'code_length', registry=ds.unit_registry)

tng.ray_func(ds, start_position=start, end_position=end, sn_ratio=18, data_filename=ray_filename+'.h5', complete_filename=f'{ray_filename}.hdf5', spectral_filename=ray_filename+'_spectra')
		
#loaded_ray_h5=h5py.File('test_file_complete.hdf5', 'r+')
#loaded_ray=yt.load('test_file_complete.hdf5')

#impact_param=loaded_ray._ray('impact_param')
#g130mflux=loaded_ray._ray('flux', 'COS-G130M')