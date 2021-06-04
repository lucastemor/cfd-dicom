'''
Read xdmf, voxelize based on resolution
Important note that slicing is no longer arbitrary, so need to read land slice along lowest coord dimension
'''

import vtk 
import pyvista as pv 
import h5py 
import numpy as np
from pathlib import Path

VOXEL_SIZE = 0.2    #in mm, assuming mesh lengths are also in mm
N_STEPS = 100       #number of steps to sample. Will take this many steps uniformly over all available steps

CLIP = True

CLIP_BOX = [19.8571, 33.8062,
			3.79448, 22.7821,
			91.6152, 107.278] #sac and parent for MCA 07

if __name__ == '__main__':
	
	case_name = 'MCA07'
	data_path = Path('/Volumes/lucas-ssd/MASc/Ubuntu_backup/dicom/pipe_ipcs_ab_cn_4D027_MCA07_constant_ts9600_cycles2_uOrder1')
	mesh_file = data_path/'4D027_MCA07.h5'

	outdir = Path('/Users/lucas/Documents/School/BSL/Horos/isotropic_voxels')

	all_timesteps = sorted(data_path.glob('*ts=*'))
	sample_indices = np.round(np.linspace(0, len(all_timesteps) - 1, N_STEPS)).astype(int)
	timesteps_to_read = [all_timesteps[i] for i in sample_indices]	

	with h5py.File(mesh_file, 'r') as hf:
		points = np.array(hf['Mesh']['coordinates'])
		cells = np.array(hf['Mesh']['topology'])

		celltypes = np.empty(cells.shape[0], dtype=np.uint8)
		celltypes[:] = vtk.VTK_TETRA

		offset = np.arange(0,cells.shape[0],cells.shape[1])   

		cell_type = np.ones((cells.shape[0], 1), dtype=int) * 4
		cells = np.concatenate([cell_type, cells], axis = 1)
		mesh = pv.UnstructuredGrid(cells.ravel(), celltypes, points)
		surf = mesh.extract_surface()

		mesh.point_arrays['ids'] = np.arange(0,mesh.n_points,1)

	if CLIP:
		mesh = mesh.clip_box(bounds=CLIP_BOX, invert=False)

	bounding_box = mesh.bounds
	x_length = np.abs(bounding_box[0] - bounding_box[1])
	y_length = np.abs(bounding_box[2] - bounding_box[3])
	z_length = np.abs(bounding_box[4] - bounding_box[5])
	
	n_voxels_x = int(np.floor(x_length/VOXEL_SIZE))
	n_voxels_y = int(np.floor(y_length/VOXEL_SIZE))
	n_voxels_z = int(np.floor(z_length/VOXEL_SIZE))
	
	for i,timestep_file in enumerate(timesteps_to_read):
		print(f'Step {i} of {len(timesteps_to_read)}',end = "\r")
		hf = h5py.File(timestep_file, 'r')
		u = np.array(hf['Solution']['u'])
		mesh.point_arrays['u'] = u[mesh.point_arrays['ids']]

		# Compute Q
		mesh = mesh.compute_derivative(qcriterion=True,gradient=False,scalars='u')

		structured_mesh = pv.create_grid(mesh, dimensions=(n_voxels_x,n_voxels_y,n_voxels_z))
		resampled_image = structured_mesh.sample(mesh)
		resampled_image.save(outdir/case_name/f"resampledimage_{str(i).zfill(4)}.vti")
