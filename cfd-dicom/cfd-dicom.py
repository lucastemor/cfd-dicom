'''
Command line tool fo converting unstructured CFD data to DICOM series


Please see README.md for description of args

'''


import sys,os, re, h5py
import numpy as np

import pyvista as pv

from vtk import VTK_TETRA
from pathlib import Path

import write_dicom as wd
import preprocessing

from scipy.ndimage import gaussian_filter

def get_scaling_function(quantity):
	if quantity == 'u':
		scaling_function = preprocessing.get_umag_pixel_array
	elif quantity == 'q':
		scaling_function = preprocessing.get_q_pixel_array
	else:
		print (f'quantity {quantity} unknown or unspecified')
		sys.exit()	

	return scaling_function


if __name__ == '__main__':

	########### Read args ##############
	
	mesh_path  = Path('/Users/BSL/Documents/AneurysmData/c0053_ACA_T/c0053_ACA_T_mesh.h5') #sys.argv[1]  
	n_steps    = 1200 #sys.argv[2]  
	voxel_size = 0.1 #sys.argv[3] 
	quantities = ['path'] #sys.argv[4:]

	########## For pathlines #########
	particle_polydata_path = Path('/Users/BSL/Documents/AneurysmData/ComputedPathlines/c0053/500randomseed')
	sigma 		 = 0.1
	track_length = 1
	step_stride  = 2
	particle_stride
	####################################

	######## Prepare spatiotemporal sampling, output, and mesh clipping #############

	case_name = mesh_path.parent.name + f'_voxsize={voxel_size}_nsteps={n_steps}'
	outdir = Path('./output').joinpath(case_name)

	if os.path.exists(outdir) == False:	
		os.mkdir(outdir)

	all_timesteps = sorted(mesh_path.parent.glob('*ts=*'))
	sample_indices = np.round(np.linspace(0, len(all_timesteps) - 1, n_steps)).astype(int)
	timesteps_to_read = [all_timesteps[i] for i in sample_indices]	

	with h5py.File(mesh_path, 'r') as hf:
		points = np.array(hf['Mesh']['coordinates'])
		cells = np.array(hf['Mesh']['topology'])

		celltypes = np.empty(cells.shape[0], dtype=np.uint8)
		celltypes[:] = VTK_TETRA

		offset = np.arange(0,cells.shape[0],cells.shape[1])   

		cell_type = np.ones((cells.shape[0], 1), dtype=int) * 4
		cells = np.concatenate([cell_type, cells], axis = 1)
		mesh = pv.UnstructuredGrid(cells.ravel(), celltypes, points)
		surf = mesh.extract_surface()

		mesh.point_data['ids'] = np.arange(0,mesh.n_points,1)

	p = pv.Plotter()
	p.add_mesh_clip_box(surf, color='white')
	p.add_text('Select geometry ROI. Close window when finished', position='upper_left', font_size=18)
	p.show()	

	bounding_box = p.box_clipped_meshes[0].extract_surface().bounds
	mesh = mesh.clip_box(bounds=bounding_box, invert=False)

	x_length = np.abs(bounding_box[0] - bounding_box[1])
	y_length = np.abs(bounding_box[2] - bounding_box[3])
	z_length = np.abs(bounding_box[4] - bounding_box[5])
	
	n_voxels_x = int(np.floor(x_length/voxel_size))
	n_voxels_y = int(np.floor(y_length/voxel_size))
	n_voxels_z = int(np.floor(z_length/voxel_size))

	##################################################################################



	#################### DICOM stack initialize ##########################

	dicom_series = {}
	for quantity in quantities:
		series_name = f'{quantity}_' + case_name
		series_path = outdir.joinpath(quantity)

		if os.path.exists(series_path) == False:	
			os.mkdir(series_path)

		stack = wd.dicom_stack(write_dir=series_path,n_timesteps=len(timesteps_to_read),case_name=series_name)	
		stack.quantity = quantity
		dicom_series[quantity] = stack

	#################################################################################

	
	################# For scalars ############################
	for i,timestep_file in enumerate(timesteps_to_read):
		print(f'Scalars: step {i} of {len(timesteps_to_read)-1}',end = "\r")
		
		if 'u' in quantities or 'q' in quantities:
			hf = h5py.File(timestep_file, 'r')
			u = np.array(hf['Solution']['u'])
			mesh.point_data['u'] = u[mesh.point_data['ids']]

			if 'q' in quantities:
				mesh = mesh.compute_derivative(qcriterion=True,gradient=False,scalars='u')

			structured_mesh = pv.create_grid(mesh, dimensions=(n_voxels_x,n_voxels_y,n_voxels_z))
			resampled_image = structured_mesh.sample(mesh)

		for quantity in quantities:
			if quantity != 'path':
				scaling_function = get_scaling_function(quantity)
				array_3d = scaling_function(resampled_image)
				dicom_series[quantity].write_isotemporal_slices(array_3d)
				#resampled_image.save(outdir/f"resampledimage_{str(i).zfill(4)}.vti")
	##################################################################################


	#################### Pathlines ##########################
	if 'path' in quantities:
		particle_time_series_files = sorted(list(particle_polydata_path.glob('*.vtp')),key = lambda x : int(re.findall(r'([\d]+).{}'.format('vtp'),str(x))[0]))[::step_stride]
		photographed_frames = np.arange(track_length,len(particle_time_series_files),track_length)
		
		if len(photographed_frames) > n_steps:
			idx = np.round(np.linspace(0, len(photographed_frames) - 1, n_steps)).astype(int)
			frames_to_write = photographed_frames[idx]
		else:
			frames_to_write = photographed_frames
		
		dicom_series['path'].n_timesteps = len(frames_to_write)

		#have to read in an arbitrary sample mesh to get geo
		hf = h5py.File(timesteps_to_read[0], 'r')
		u = np.array(hf['Solution']['u'])
		mesh.point_data['u'] = u[mesh.point_data['ids']]
		mesh.point_data['u'] = u[mesh.point_data['ids']]
		structured_mesh = pv.create_grid(mesh, dimensions=(n_voxels_x,n_voxels_y,n_voxels_z))
		resampled_image = structured_mesh.sample(mesh)
		scaling_function = get_scaling_function('u')
		array_3d = scaling_function(resampled_image)

		model_overlay = np.where(array_3d>0,1,0)

		voxel_array = np.zeros_like(model_overlay)	

		voxel_size = np.max(resampled_image.spacing)
		sigma = sigma/voxel_size
		
		track_counter = 1
		for t,timestep in enumerate(particle_time_series_files[:-1]):
			print (f'\n Pathlines: step {t} of {len(particle_time_series_files)-1}', end = '\r')
			track_counter +=1

			next_step = particle_time_series_files[t+1]

			all_particles = pv.read(timestep)
			next_all_particles = pv.read(next_step)

			common_ids = np.intersect1d(all_particles.point_data['ParticleId'],next_all_particles.point_data['ParticleId'])

			for particle_id in common_ids:

				point = all_particles.points[np.where(all_particles.point_data['ParticleId'] == particle_id)][0]
				next_point = next_all_particles.points[np.where(next_all_particles.point_data['ParticleId'] == particle_id)][0]

				#point = particle.points[0]
				#next_point = next_particle.points[0]

				#find closest point on grid
				point_id_on_grid = resampled_image.FindPoint(point)
				next_point_id_on_grid = resampled_image.FindPoint(next_point)

				if point_id_on_grid != -1 and next_point_id_on_grid !=-1:

					#we have a valid point

					x1,y1,z1 = preprocessing.point_id_to_structured_coords(resampled_image, point_id_on_grid)

					x2,y2,z2 = preprocessing.point_id_to_structured_coords(resampled_image, next_point_id_on_grid)
					
					voxel_spline = preprocessing.Bresenham3D(x1,y1,z1,x2,y2,z2)

					for line_point in voxel_spline:
						z=line_point[0] #flip for writing from np array to dcm directly
						y=line_point[1]
						x=line_point[2]

						voxel_array[x,y,z] +=1
				
				else:
					print (f'skipped particle {particle_id} at step {t} because not in domain',end="\r")
				

			if t in photographed_frames and t in frames_to_write:
				blurred = gaussian_filter(voxel_array,sigma=sigma)
				voxel_array = blurred/blurred.max()

				'''
				#see https://vtk.org/Wiki/VTK/Writing_VTK_files_using_python
				#dimensions
				nx, ny, nz = vti.GetDimensions()
				ncells = nx * ny * nz
				npoints = (nx + 1) * (ny + 1) * (nz + 1) 

				outfile = f'/mnt/3414B51914B4DED4/dicom/data/voxelize_pathlines/500randomseed/vti/TEST3dbresenhm_200x200x200_track{track_length}_blurred_sigma{sigma}_{str(t).zfill(4)}'

				imageToVTK(outfile, origin= vti.GetOrigin(), spacing= vti.GetSpacing(), pointData = {"pathline" : blurred} ) #can work for cell data too
				'''
				
				
				voxel_array *= 1000
				voxel_array+=model_overlay
				#voxel_array = np.where(voxel_array==0,model_overlay,voxel_array)

				#print (f'-------------Step {t}-------------------')
				#print (voxel_array.min(),voxel_array.max())
				#print (voxel_array.astype(np.uint16).min(),voxel_array.astype(np.uint16).max())
				dicom_series['path'].write_isotemporal_slices(voxel_array)
				

				#reset voxel to 0 at every timstep. Or in future use 4d arrray and write out slices?
				voxel_array = np.zeros_like(model_overlay)	
				track_counter = 1

			elif t in photographed_frames:
				voxel_array = np.zeros_like(model_overlay)	
				track_counter = 1


			#print (f'Failed at time {t}')
		


			#new_point = pv.PolyData(np.array(vti.GetPoint(vti.FindPoint(point))))
			#new_point.field_arrays['TimeValue'] =  particle.field_arrays['TimeValue']


			#new_point.save(f'/mnt/3414B51914B4DED4/dicom/data/voxelize_pathlines/SINGLE_POINT/single_point_tracked_polydata/200x200grid_new_point/resampled_point_{str(t).zfill(4)}.vtp')
