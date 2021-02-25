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

from pyevtk.hl import imageToVTK

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
	vti_series_path = '/home/lucas/Documents/viz/renders/Horos/surgical/MCA07/velocity_200x200/'
	png_series_path = '/Users/lucas/Documents/School/BSL/cfd-dicom/foo/pathline/'

	particle_polydata_path = '/mnt/3414B51914B4DED4/dicom/data/voxelize_pathlines/SINGLE_POINT/single_point_tracked_polydata/'


	vti_time_series_files = glob_sort_files(vti_series_path,'vti')
	particle_time_series_files = glob_sort_files(particle_polydata_path,'vtp')

	reader = vtk.vtkXMLImageDataReader()
	reader.SetFileName(vti_time_series_files[0])
	reader.Update()
	vti = reader.GetOutput()

	voxel_array = np.zeros(vti.GetDimensions())

	#outdir = '/mnt/3414B51914B4DED4/dicom/data/voxelize_pathlines/SINGLE_POINT/single_point_gaussian_dicom_test/'
	#dicom_stack = wd.dicom_stack(write_dir=outdir,n_timesteps=len(time_series_files),case_name=case_name)

	sigma=1
	track_length=5  #similar to shutter speed, or better track length from temporal particles to pathlines in paraview

	track_counter = 1
	for t,timestep in enumerate(particle_time_series_files[:-1]):

		next_step = particle_time_series_files[t+1]

		particle = pv.read(timestep)
		next_particle = pv.read(next_step)

		if len(next_particle.points>0):

			point = particle.points[0]
			next_point = next_particle.points[0]

			#find closest point on grid
			point_id_on_grid = vti.FindPoint(point)
			if point_id_on_grid != -1:

				#we have a vild point
				track_counter +=1

				x1,y1,z1 = point_id_to_structured_coords(vti, point_id_on_grid)

				next_point_id_on_grid = vti.FindPoint(next_point)
				x2,y2,z2 = point_id_to_structured_coords(vti, next_point_id_on_grid)
				
				voxel_spline = Bresenham3D(x1,y1,z1,x2,y2,z2)

				for line_point in voxel_spline:
					x = line_point[0]
					y=line_point[1]
					z=line_point[2]

					voxel_array[x,y,z] = 1
				
				if track_counter == track_length:
					blurred = gaussian_filter(voxel_array,sigma=sigma)

					#see https://vtk.org/Wiki/VTK/Writing_VTK_files_using_python
					#dimensions
					nx, ny, nz = vti.GetDimensions()
					ncells = nx * ny * nz
					npoints = (nx + 1) * (ny + 1) * (nz + 1) 

					outfile = f'/mnt/3414B51914B4DED4/dicom/data/voxelize_pathlines/SINGLE_POINT/single_point_track1_vti_test/3dbresenhm_200x200x200_track{track_length}_blurred_sigma{sigma}_{str(t).zfill(4)}'

					imageToVTK(outfile, origin= vti.GetOrigin(), spacing= vti.GetSpacing(), pointData = {"pathline" : blurred} ) #can work for cell data too

					#reset voxel to 0 at every timstep. Or in future use 4d arrray and write out slices?
					voxel_array = np.zeros(vti.GetDimensions())

					track_counter = 1

			
			else:
				print (f'skipped particle at step {t} because not in domain')

		else:
			print (f'skipped timestep {t} becasue no particle')


		#print (f'Failed at time {t}')
	


		#new_point = pv.PolyData(np.array(vti.GetPoint(vti.FindPoint(point))))
		#new_point.field_arrays['TimeValue'] =  particle.field_arrays['TimeValue']


		#new_point.save(f'/mnt/3414B51914B4DED4/dicom/data/voxelize_pathlines/SINGLE_POINT/single_point_tracked_polydata/200x200grid_new_point/resampled_point_{str(t).zfill(4)}.vtp')

	