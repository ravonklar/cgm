#For now, most notes on local data are kept in tng_localdata_tutorial.py, which isn't necessarily meant to run as an actual file. But below are some important notes:

#IllustrisTNG file preamble:

#Particle types:

#0) gas
#1) dm
#2) not used
#3) tracer
#4) stars
#5) bh


#SIMULATION UNITS:

#Coordinates: ckpc/h, comoving--> scale calculated later converts
#Mass: 10^10/h solar masses
#Density: 10^10/h solar masses/comoving (kpc/h^3)
#Internal energy (per unit mass): (km/s)^2
#Velocity: km/s*sqrt(a)
#NOTE: some of these parameters are scaled by h, others are not!
#NOTE ON UNITS IN GENERAL: It seems that if you have distance in kpc/h, and you want it in kpc, you should multiply by h.
#But! Consider: if you have mass in units of 10^10Msun, and you have m=4, then m=4*(10^10Msun). Similarly, if units are 10^10Msun/h, and m=4, m=4*10^10Msun/h! Ie, if distance is 5, and units are ckpc/h, then distance=5*ckpc/h.
#MW parameters

#MW total mass: ~0.8-1.5e12=80-150*10**10
#MW SFR: 0.68-1.45 M_sun/year

#COMOVING VS PROPER DISTANCE:
#The simulation units are in comoving coordinates and km*sqrt(a)/s (so we multiply by sqrt(a) to get km/s). These do not change as the Universe expands. We have been translating this into *proper* distance and velocity, which do change with the expansion. But they are the true physical distances.
#Right now I don't think the distinction matters much, because we are looking at interacts between particles within a galaxy-- their relative positions/velocities should not be altered by Hubble flow.


#"CHARMING" QUIRKS: Dealing with datasets with one value can be annoying. Sometimes the value of, say, subinfo['phi'] is returned as "<HDF5 dataset "vairable_name": shape (), type "<f8">", while subinfo['phi'][0] returns an error about illegal slicing. To see the actual value do

#subinfo['phi'][()]    (subinfo['phi'].value also works, but is deprecated)
#A scalar value like one angle isn't physically stored in memory the same way a list is (I think), so you can't just have subinfo['keyword'] return an int or float, it has to be an array-like object to be in an HDF5 format. But it *is* just one value, and that can be physically meaningful to us as scientists (ie, you may have a list of x-coords for stars. This is normally an array, but you just so happen to have only one star in list. This should be stored differently from metadata about the subhalo, like stellar half-mass, where you expect 1 number)
#ALSO! IMPORTANT: THE HIGHEST LEVEL OF DATA IN HDF5 FILE IS NOT THE h! To see the highets level of data, do f=h5py.File()-->f.keys(), Header is simply one of the keys. No idea what all the others do. 
import numpy as np
import yt
yt.enable_plugins()
import copy
from unyt import unyt_array
import trident
import h5py
import six
import requests
import pickle
import astropy
#from illustris_python.util import partTypeNum
from os.path import isfile
from astropy import constants as c, units as u
from astropy.cosmology import Planck15 as cosmo, z_at_value


#slis 230M, COS: G185M, G225M, G285M, G230L
#can see various version of redshift: redshift=cosmo, redshift_dopp=peculiar, redshift_eff=total
cos_grisms=['COS-G130M', 'COS-G160M', 'COS-G140L', 'COS-G185M', 'COS-G225M', 'COS-G285M', 'COS-G230L']
grism_info=[[899, 1469, 0.00997, 'avg_COS_G130M.txt'], [1342, 1798, 0.01223, 'avg_COS_G160M.txt'], [815, 2391, 0.0803, 'avg_COS_G140L.txt'], [1670, 2127, 0.037, None], [2070, 2527, 0.033, None], [2480, 3229, 0.04, None], [1334, 2361, 0.39, None]]
instruments={cos_grisms[i]:grism_info[i] for i in range(len(cos_grisms))}
#sgs=[trident.SpectrumGenerator(lambda_min=grism_info[i][0], lambda_max=grism_info[i][1], dlambda=grism_info[i][2], line_database='lines.txt') for i in range(3)]


ds_index_holder, ray_index_holder, common_values_holder=0, 0, 0


basepath='./TNG100-1/output'
h=0.6774
H0=100*h*u.km.to('Mpc')
G=c.G.to('kpc3 / (solMass s2)').value*10**10  #this converts G to correspond with Illustris units (aside from a factor of h) 

def temp_calc(u, x_e):
  X_H=0.76 #"Hydrogen mass function"
  k_B=1.380648e-16 #in CGS (ergs)
  UEUM=(3.086e21/3.1536e16)**2 #unit energy/unit mass. It is equal to (unitLength/unitTime)**2=(kpc/Gyr)**2, in CGS. 3.086 cm in kpc, 3.1536e16s in Gyr 
  mp=1.6726e-24 #proton mass in grams
  gamma=5.0/3.0
  mu=4./(1+3*X_H+4*X_H*x_e)*mp #mean molecular weight IN GRAMS
  T=(gamma-1.)*(u/k_B)*UEUM*mu
  return T, mu


groups=['Config', 'Header', 'halo_and_sub_properties', 'ray_properties', 'PartType5', 'PartType0']
ray_prop_names=['ray_uv', 'impact_uv', 'impact_parameter']
def ray_func(ds, start_position, end_position, instruments=instruments, sn_ratio=18,
             lines='all', ftype="PartType0", fields=None,
             solution_filename=None, data_filename=None, spectral_filename=None, complete_filename=None,
             trajectory=None, redshift=None, field_parameters=None,
             setup_function=None, load_kwargs=None,
             line_database='MortCash0.txt', ionization_table=None, interactive=False):
           
  origin=unyt_array(75000/2, units='code_length', registry=ds.unit_registry)
  #data_filename is for the ray.h5 object; spectral filename is for the images/data files; complete_filename is for the final hdf5 file we make here     
  #start and end positions are assumed to be input as either unyt_arrays or bare lists/arrays. If bare, assume in simulation coords. First thing we do is unscale the start and end points from kpc to code lengths

  if hasattr(start_position, 'units'):
    start_position, end_position=start_position.to('code_length'), end_position.to('code_length')
  else:
    start_position, end_position=unyt_array(start_position, units='code_length', registry=ds.unit_registry), unyt_array(end_position, units='code_length', registry=ds.unit_registry)
    
    #Now we find the unit vector of our ray, the impact parameter, and the unit vector from the origin to the impact point. We do all this before shifting the ray position to fit where yt thinks the actual particles are
  ray_vec=start_position-end_position
  ray_uvec=ray_vec.value/np.linalg.norm(ray_vec)  #unit vectors are ironically unitless
  distance_along_ray=np.sum(-start_position*ray_uvec)  #dot product tells us how far along ray the point closest to the origin is
  x=start_position+distance_along_ray*ray_uvec  #point along ray closest to origin; because origin=galactic center, can do next line easy:
  ip=unyt_array(np.linalg.norm(x), units='code_length', registry=ds.unit_registry)   #this is ip in code length units, but comes out of norm() unitless
  uv_to_ip=x.value/ip.value
    
  print('print statement 1')
  #Now we make the actual ray
  ray=trident.make_simple_ray(ds, start_position=start_position+origin, end_position=end_position+origin, data_filename=data_filename, lines=lines, ftype='PartType0', line_database=line_database, fields=[('PartType0', 'ParticleIDs')])
  dl=copy.deepcopy(ray.r['gas', 'dl'])
  

  
  ray_props=[ray_uvec, uv_to_ip, ip]
  
  
  #Everything has been calculated; now we just need to organize/store data in hdf5 file    
  #First we copy/store metadata in a new file
  with h5py.File(complete_filename, 'w') as f:
    for grp in groups:
      f.create_group(grp)
    f['Config'].attrs['VORONOI'], f['Config'].attrs['RAY']=1, 1  
    
  #trying to move away from using attrs, but this set-up is essential to reloading the data files, as yt assumes this arrangement. Same with header below
    for key, value in ds.headvals.items():
      if key=='NumPart_ThisFile':
        f['Header'].attrs[key]=np.int32([len(ray.r['gas', 'ParticleIDs']), 0., 0., 0., 0., ds.bh('count')])
      else:
        f['Header'].attrs[key]=value
    for key, value in ds.hsvals.items():
      f['halo_and_sub_properties'].create_dataset(key, data=value)
  
  #Now we store data about the ray itself 
    for i in range(len(ray_props)):
      f['ray_properties'].create_dataset(ray_prop_names[i], data=np.array(ray_props[i]))
    #f['ray_properties'].create_dataset('instrument_order', data=np.array(list(instruments.keys())))  #tells us the order of the instruments
    #We use our instruments to create spectral data from the ray
    lambda_, tau, flux, error=[], [], [], []
    for inst_name, inst_props in instruments.items():   #generate/save spectra for different instruments
      sg=trident.SpectrumGenerator(lambda_min=inst_props[0], lambda_max=inst_props[1], dlambda=inst_props[2], line_database=line_database)
      sg.make_spectrum(ray, lines=lines)
      if isinstance(inst_props[3], str):
        sg.apply_lsf(filename=inst_props[3])
      else:
        sg.apply_lsf(function='gaussian', width=4)     #until we get actual LSF files
      sg.add_gaussian_noise(sn_ratio)
      f['ray_properties'].create_dataset(inst_name, data=np.array([sg.lambda_field, sg.tau_field, sg.flux_field, sg.error_field]))
      if interactive==True:
        sg.save_spectrum(spectral_filename+'_'+inst_name+'.txt')
        sg.plot_spectrum(title=spectral_filename+' '+inst_name, filename=spectral_filename+'_'+inst_name+'.png')
      
      
    f['ray_properties'].create_dataset('property_order', data=np.array([b'lambda', b'tau', b'relative_flux', b'flux_error']))
    #f['ray_properties']['lambda'].create_dataset(key, data=np.array(lambda_))
    #f['ray_properties']['tau'].create_dataset(key, data=np.array(tau))
    #f['ray_properties']['flux'].create_dataset(key, data=np.array(flux))
  
   
  #Now we store particle data. First, create mask for gas data in ds object
    common_values, ds_index, ray_index=np.intersect1d(ds.gas('ParticleIDs'), ray.r['gas', 'ParticleIDs'], assume_unique=True, return_indices=True)
    l_ray=ray.r['gas', 'l'][ray_index]  #this is a list of l values in the same order as ray properties will be in. We want to make a list of indeces for the ray and for the ds that go in order of increasing 'l'. We can zip these values with ds_index and ray_index to get a new order for both which goes in order of 'l' instead of randomly
    placeholder1=list(zip(l_ray, ds_index, ray_index))
    placeholder1.sort()  #sorts indeces by corresponding 'l' value
    placeholder2=list(zip(*placeholder1))
    common_values, ds_index, ray_index=np.array(placeholder2[0]), np.array(placeholder2[1]), np.array(placeholder2[2])  #These are now the right order of indeces to give particles in order along the ray
    global ds_index_holder, ray_index_holder, common_values_holder
    ds_index_holder, ray_index_holder, common_values_holder=ds_index, ray_index, common_values
    #To keep all data in right order, need to use these "masks" on both ray object and ds object. Even though all ray object particles are used, they aren't in the same order (ie, ray_index is not just [0, 1, 2, 3....]). 
    for field in ds.derived_field_list:
      if 'PartType5' in field:  #if field tuple starts with PartType5
        f['PartType5'].create_dataset(field[1], data=ds.bh(field[1]))
      elif 'PartType0' in field: 
        if 'count' in field:
          f['PartType0'].create_dataset(field[1], data=len(ray_index))
        else:
          print(field)
          f['PartType0'].create_dataset(field[1], data=ds.gas(field[1])[ds_index])
    #We take properties from the dataset first; some auto-generated properties in the ray may have the same name as properties I made up in the dataset. My own properties are better, so we take them. If we find the same property in the ray we just skip it
    for field in ray.derived_field_list:
      if field==('gas', 'dl'):
        f['PartType0'].create_dataset(field[1], data=dl[ray_index])
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
      elif 'gas' in field and (field[1]=='l' or 'number_density' in field[1] or 'ion_fraction'in field[1] or 'nuclei_density' in field[1]):
        try:
          f['PartType0'].create_dataset(field[1], data=ray.r[field][ray_index])
        except OSError:   #this error is thrown if we try to create a dataset with a name that is already used. Should only happen if we grab a property from ray that already existed for all particles in ds
          print(field, 'OSError exception')
        except:
          print(field, 'This should only happen with hydrogen nuclei density and neutral hydrogen density')
          
  
  #Finally done storing/organizing data; now just save hdf5 file. Note that we are saving several properties from ds that are normally calculated from TNG properties on disk when loaded with yt; for ray files, these properties will already be calculated and saved to disk. When you use yt.load() on them the properties are needlessly recalculated, but nothing bad happens
    f.close()
  if interactive==True:
    return yt.load(complete_filename)
  else:
    return complete_filename
    
        
    
               
            
            
  

def ray_transform(theta, phi, origin, ray): #ray object should be passed as [ray_start, ray_end]
  ray_start, ray_end=np.array(ray[0]).reshape(3, 1), np.array(ray[1]).reshape(3, 1)
  R_zT=np.array([[np.cos(phi), np.sin(phi), 0], [-1*np.sin(phi), np.cos(phi), 0], [0, 0, 1]])
  R_yT=np.array([[np.cos(theta), 0, -1*np.sin(theta)], [0, 1, 0], [np.sin(theta), 0, np.cos(theta)]])
  rotated_start, rotated_end=(R_zT @ R_yT @ ray_start).reshape(1, 3), (R_zT @ R_yT @ ray_end).reshape(1, 3)
  return np.array([rotated_start[0]+origin, rotated_end[0]+origin])




def gasorient(gas, origin, groupv, radial_cut):
#this function takes the velocities, masses, etc of gas particles around and origin and finds the orientation of the cool gaseous disk, if any, by using angular momentum. Ideally, pass the gas+origin of central galaxy in halo only. It needs a radial cut, currently using stellar half-mass radius
  if not isinstance(origin, np.ndarray):
    origin=np.array(origin, dtype=np.float32)
  if not isinstance(groupv, np.ndarray):
    groupv=np.array(groupv, dtype=np.float32)
  T_cut=1e4
  dxdydz=gas['Coordinates']-origin #dxdydz[:, 0]=dx, etc
  vxvyvz=gas['Velocities']-groupv
  r2=np.sum((dxdydz)**2, axis=1) #NOTE r IS SQUARED STILL. Multiply by (scale/h)**2 for physical units. But because we take trig ratios, can use code units throughout. 
  T, mu=temp_calc(gas['InternalEnergy'][:], gas['ElectronAbundance'][:])
  mask=(r2<radial_cut**2) & (T<T_cut)
#Now we calculate angular momentum of central cool gas about each axis
  #dxdydz, vxvyvz, mass=dxdydz[mask], vxvyvz[mask], gas['Masses'][mask]
  Lx=np.sum(gas['Masses'][mask]*(dxdydz[mask, 1]*vxvyvz[mask, 2]-dxdydz[mask, 2]*vxvyvz[mask, 1]))
  Ly=np.sum(gas['Masses'][mask]*(dxdydz[mask, 2]*vxvyvz[mask, 0]-dxdydz[mask, 0]*vxvyvz[mask, 2]))
  Lz=np.sum(gas['Masses'][mask]*(dxdydz[mask, 0]*vxvyvz[mask, 1]-dxdydz[mask, 1]*vxvyvz[mask, 0]))
#Finally, the orientation:
  rot_theta, rot_phi=np.arctan2(np.sqrt(Lx**2+Ly**2), Lz), np.arctan2(Ly, Lx)
  return rot_theta, rot_phi
  
  




      
    



























