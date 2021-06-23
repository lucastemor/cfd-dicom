import vtk,sys
import numpy as np
import glob,re,os

import preprocessing
import write_dicom as wd

from skimage import io 
from skimage.color import rgb2gray 

import matplotlib.pyplot as plt
import pyvista as pv

from scipy.ndimage import gaussian_filter

#from pyevtk.hl import imageToVTK

def glob_sort_files(series_path,extension):
	return sorted(glob.glob(series_path+f'/*.{extension}'),key = lambda x : int(re.findall(r'([\d]+).{}'.format(extension),x)[0]))#[::5]

def Bresenham3D(x1, y1, z1, x2, y2, z2): 
	"""
	# Python3 code for generating points on a 3-D line  
	# using Bresenham's Algorithm 
	Assumes unit grid spacing (i.e., structured coords)
	Source https://www.geeksforgeeks.org/bresenhams-algorithm-for-3-d-line-drawing/
	"""
	ListOfPoints = [] 
	ListOfPoints.append((x1, y1, z1)) 
	dx = abs(x2 - x1) 
	dy = abs(y2 - y1) 
	dz = abs(z2 - z1) 
	if (x2 > x1): 
		xs = 1
	else: 
		xs = -1
	if (y2 > y1): 
		ys = 1
	else: 
		ys = -1
	if (z2 > z1): 
		zs = 1
	else: 
		zs = -1

	# Driving axis is X-axis" 
	if (dx >= dy and dx >= dz):         
		p1 = 2 * dy - dx 
		p2 = 2 * dz - dx 
		while (x1 != x2): 
			x1 += xs 
			if (p1 >= 0): 
				y1 += ys 
				p1 -= 2 * dx 
			if (p2 >= 0): 
				z1 += zs 
				p2 -= 2 * dx 
			p1 += 2 * dy 
			p2 += 2 * dz 
			ListOfPoints.append((x1, y1, z1)) 

	# Driving axis is Y-axis" 
	elif (dy >= dx and dy >= dz):        
		p1 = 2 * dx - dy 
		p2 = 2 * dz - dy 
		while (y1 != y2): 
			y1 += ys 
			if (p1 >= 0): 
				x1 += xs 
				p1 -= 2 * dy 
			if (p2 >= 0): 
				z1 += zs 
				p2 -= 2 * dy 
			p1 += 2 * dx 
			p2 += 2 * dz 
			ListOfPoints.append((x1, y1, z1)) 

	# Driving axis is Z-axis" 
	else:         
		p1 = 2 * dy - dz 
		p2 = 2 * dx - dz 
		while (z1 != z2): 
			z1 += zs 
			if (p1 >= 0): 
				y1 += ys 
				p1 -= 2 * dz 
			if (p2 >= 0): 
				x1 += xs 
				p2 -= 2 * dz 
			p1 += 2 * dy 
			p2 += 2 * dx 
			ListOfPoints.append((x1, y1, z1)) 
	return ListOfPoints 



def point_id_to_structured_coords(vti, point_id):
	"""
	vti written in paraview will count from (0,0,0) to (dim,dim,dim) incrementing xdime+1 every time and periodically loopiung back to x=0 but then iuncrmeenting y up to y dim .. etc..
	"""

	x_dim,y_dim,z_dim = vti.GetDimensions()

	x = point_id%x_dim
	y = int(((point_id-x)/x_dim)%y_dim)
	z = int((point_id-x-y*y_dim)/(x_dim*y_dim))

	return x,y,z

if __name__=='__main__':
	"""
	Series id tracks "global" posistion in the dicom stack
	Slice id tracks "local" position in the timestep stack
	"""
	n_DICOM_frames = 100 #max in horos is 100 

	track_length= 4  #similar to shutter speed, or better track length from temporal particles to pathlines in paraview
	step_stride = 1
	sigmas = [0,0.05,0.1]
	voxel_sizes = [0.2,0.1,0.05]
	#densities = [500,250,125,63]

	vti_series_path = '/Users/lucas/Documents/School/BSL/Horos/isotropic_voxels/comparison_voxelsize/'
	particle_polydata_path = '/Volumes/lucas-ssd/MASc/Ubuntu_backup/dicom/200mask/polydata/'

	file_paths = glob.glob(vti_series_path+'*.vti')

	for voxel_path in file_paths:

		reader = vtk.vtkXMLImageDataReader()
		reader.SetFileName(voxel_path)
		reader.Update()
		vti = reader.GetOutput()

		for sigma in sigmas:
			case_name = 'SEED_compare_'+os.path.splitext(os.path.split(voxel_path)[1])[0]+f'_{sigma}sigma'

			outdir = f'./output/{case_name}/'
			if os.path.exists(outdir) == False:	
				os.mkdir(outdir)


			vti_time_series_files = glob_sort_files(vti_series_path,'vti')
			particle_time_series_files = glob_sort_files(particle_polydata_path,'vtp')[::step_stride]

			photographed_frames = np.arange(track_length,len(particle_time_series_files),track_length)

			if len(photographed_frames) > n_DICOM_frames:
				idx = np.round(np.linspace(0, len(photographed_frames) - 1, n_DICOM_frames)).astype(int)
				frames_to_write = photographed_frames[idx]
			else:
				frames_to_write = photographed_frames


			u_voxels = preprocessing.get_umag_pixel_array(vti)
			model_overlay = np.where(u_voxels>0,1,0)

			z_dim,y_dim,x_dim = vti.GetDimensions()
			voxel_array = np.zeros((x_dim,y_dim,z_dim))	

			voxel_size = np.max(vti.GetSpacing())
			sigma = sigma/voxel_size

			n_timesteps = n_DICOM_frames
			#n_timesteps= int(len(particle_time_series_files)/track_length)
			dicom_stack = wd.dicom_stack(write_dir=outdir,n_timesteps=n_timesteps,case_name=case_name)


			track_counter = 1
			for t,timestep in enumerate(particle_time_series_files[:-1]):

				track_counter +=1

				next_step = particle_time_series_files[t+1]

				all_particles = pv.read(timestep)
				next_all_particles = pv.read(next_step)

				common_ids = np.intersect1d(all_particles.point_arrays['ParticleId'],next_all_particles.point_arrays['ParticleId'])

				for particle_id in common_ids:

					point = all_particles.points[np.where(all_particles.point_arrays['ParticleId'] == particle_id)][0]
					next_point = next_all_particles.points[np.where(next_all_particles.point_arrays['ParticleId'] == particle_id)][0]

					#point = particle.points[0]
					#next_point = next_particle.points[0]

					#find closest point on grid
					point_id_on_grid = vti.FindPoint(point)
					next_point_id_on_grid = vti.FindPoint(next_point)

					if point_id_on_grid != -1 and next_point_id_on_grid !=-1:

						#we have a vild point

						x1,y1,z1 = point_id_to_structured_coords(vti, point_id_on_grid)

						x2,y2,z2 = point_id_to_structured_coords(vti, next_point_id_on_grid)
						
						voxel_spline = Bresenham3D(x1,y1,z1,x2,y2,z2)

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

					print (f'-------------Step {t}-------------------')
					print (voxel_array.min(),voxel_array.max())
					print (voxel_array.astype(np.uint16).min(),voxel_array.astype(np.uint16).max())
					dicom_stack.write_isotemporal_slices(voxel_array)
					

					#reset voxel to 0 at every timstep. Or in future use 4d arrray and write out slices?
					voxel_array = np.zeros((x_dim,y_dim,z_dim))	
					track_counter = 1

					break

				elif t in photographed_frames:
					voxel_array = np.zeros((x_dim,y_dim,z_dim))	
					track_counter = 1
	