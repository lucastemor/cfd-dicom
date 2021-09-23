import vtk,sys
import numpy as np
import glob,re,os

from scipy.signal import convolve2d as conv2

SCALE_FACTOR = 1000

def rescale_to_dicom(array):
	return array*SCALE_FACTOR

def struct_reshape(array_1d,dimensions_3d):
	'''
	Converting from pv structured coordinate system to numpy coordinate system -- axes are labelled differently
	'''
	z_dim, y_dim, x_dim = dimensions_3d
	return array_1d.reshape(x_dim,y_dim,z_dim)


def get_q_pixel_array(pv_structured_grid,name=None,threshold=True):
	"""
	Default name should be qcriterion from pyvista, butcan specify custom name if needed
	Can be used for other scalar quantities if cusotm name is passed
	"""

	dimensions = pv_structured_grid.dimensions

	if name==None:
		name = 'qcriterion'

	q_criterion = pv_structured_grid.point_data[name]
	if threshold: 
		q_criterion = np.where(q_criterion < 0, 0.02, q_criterion)
		q_criterion = np.where(q_criterion>3,3,q_criterion)
		#print ('Threshold!', q_criterion.min(), q_criterion.max()) 	

	q_voxel = struct_reshape(q_criterion,dimensions)

	return rescale_to_dicom(q_voxel)


def get_umag_pixel_array(pv_structured_grid,name=None):
	"""
	Array name will vary depedniong on how simulation results are saved
	here we try two possibiltiies, if not, custom name can be specifiied 
	Can be used for other vector quantities if custom name is passed
	"""

	dimensions = pv_structured_grid.dimensions

	u_mag = np.linalg.norm(np.array(pv_structured_grid.point_data['u']),axis=1)
	
	u_voxel = struct_reshape(u_mag,dimensions)

	return rescale_to_dicom(u_voxel)

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


def point_id_to_structured_coords(pv_structured, point_id):
	"""
	adaptted for pyvista sturctured grid
	vti written in paraview will count from (0,0,0) to (dim,dim,dim) incrementing xdime+1 every time and periodically loopiung back to x=0 but then iuncrmeenting y up to y dim .. etc..
	"""

	x_dim,y_dim,z_dim = pv_structured.dimensions

	x = point_id%x_dim
	y = int(((point_id-x)/x_dim)%y_dim)
	z = int((point_id-x-y*y_dim)/(x_dim*y_dim))

	return x,y,z

################### OLD functions for reading VTI ###########################


def vti_reshape(array_1d,dimensions_3d):
	'''
	Converting from vti coordinate system to numpy coordinate system -- axes are labelled differently
	'''
	z_dim, y_dim, x_dim = dimensions_3d
	return array_1d.reshape(x_dim,y_dim,z_dim)

def vti_get_q_pixel_array(vti,name=None,threshold=True):
	"""
	Default name should be qcriterion from pyvista, butcan specify custom name if needed
	Can be used for other scalar quantities if cusotm name is passed
	"""

	dimensions = vti.GetDimensions()

	if name==None:
		name = 'qcriterion'

	q_criterion = np.array(vti.GetPointData().GetArray(name))
	if threshold: 
		q_criterion = np.where(q_criterion < 0, 0.02, q_criterion)
		q_criterion = np.where(q_criterion>3,3,q_criterion)
		#print ('Threshold!', q_criterion.min(), q_criterion.max()) 	

	q_voxel = vti_reshape(q_criterion,dimensions)

	return rescale_to_dicom(q_voxel)


def vti_get_umag_pixel_array_vti(vti,name=None):
	"""
	Array name will vary depedniong on how simulation results are saved
	here we try two possibiltiies, if not, custom name can be specifiied 
	Can be used for other vector quantities if custom name is passed
	"""

	dimensions = vti.GetDimensions()

	default_names = ['u','Solution/u'] #add to this as you see fit, or can index with custom name

	if name == None:
		n_arrays = vti.GetPointData().GetNumberOfArrays()
		name = list(set([vti.GetPointData().GetArrayName(i) for i in range(n_arrays)]).intersection(default_names))[0]


	u_mag = np.linalg.norm(np.array(vti.GetPointData().GetArray(name)),axis=1)
	
	u_voxel = vti_reshape(u_mag,dimensions)

	return rescale_to_dicom(u_voxel)


'''
def psf_volume(array_3d):
	for slice_id,slice_2d in enumerate(array_3d):
		array_3d[slice_id] = psf_slice(slice_2d)
	return(array_3d)


def psf_slice(slice_2d):
	psf = np.ones((5, 5)) / 25
	blurred = conv2(slice_2d, psf, 'same')
	return(blurred)

def global_png_scale(image,global_max_pixel_value):
	return (image*(MAX_DICOM_VALUE/global_max_pixel_value) + MAX_DICOM_VALUE).astype(np.uint16)
'''



