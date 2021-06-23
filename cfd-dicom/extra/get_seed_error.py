import vtk,sys
import numpy as np
import glob,re,os

import preprocessing
import write_dicom as wd

from skimage import io 
from skimage.color import rgb2gray 

import matplotlib.pyplot as plt
import matplotlib
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

	resolutions = [25,50,75,100,150,250,350,450]
	density= 500

	errors = []

	seed_path = f'/mnt/3414B51914B4DED4/dicom/data/comparisons/seed_density_resample/seed_{density}.vtp'
	all_particles = pv.read(seed_path)
	common_ids = all_particles.point_arrays['ParticleId']

	for resolution in resolutions:
		differences = []
		voxel_path = f'/mnt/3414B51914B4DED4/dicom/data/comparisons/spatial_voxelization/velocity_{resolution}.vti'
		reader = vtk.vtkXMLImageDataReader()
		reader.SetFileName(voxel_path)
		reader.Update()
		vti = reader.GetOutput()

		for particle_id in common_ids:

			original_coords = all_particles.points[np.where(all_particles.point_arrays['ParticleId'] == particle_id)][0]

			point_id_on_grid = vti.FindPoint(original_coords)

			if point_id_on_grid != -1:

				#we have a vild point
				new_coords = vti.GetPoint(point_id_on_grid)

				differences.append(np.sqrt((original_coords[0] - new_coords[0])**2 + (original_coords[1] - new_coords[1])**2 + (original_coords[2] - new_coords[2])**2))

			else:
				print (f'skipped particle {particle_id} at step {t} because not in domain',end="\r")


		average_error = np.mean(differences)
		errors.append(average_error)



	font = {'family' : 'times','weight' : 'normal','size'   : 14}

	matplotlib.rc('font', **font)

	plt.figure(figsize=(15,10))
	plt.xlim([0,500])
	plt.ylim(0,0.5)


	plt.plot(resolutions,errors, linestyle='dashed', marker = 'o',linewidth=2,markersize=10)

	annotations = [f'({r},{np.round(e,3)})' for r,e in zip(resolutions,errors)]
	for i,label in enumerate(annotations):
		plt.text(resolutions[i],errors[i],label)

	plt.title('Mean distance error vs. voxel resolution',size=22)
	plt.xlabel('Voxel resolution (x$^3$)',size=20)
	plt.ylabel('Mean distance error [mm]',size=20)

	
	plt.savefig('/mnt/3414B51914B4DED4/dicom/paper_figs/particle_error/particle_error.png')





