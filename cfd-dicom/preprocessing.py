import vtk,sys
import numpy as np
import glob,re,os

from scipy.signal import convolve2d as conv2

SCALE_FACTOR = 1000

def rescale_to_dicom(array):
	return array*SCALE_FACTOR

def vti_reshape(array_1d,dimensions_3d):
	'''
	Converting from vti coordinate system to numpy coordinate system -- axes are labelled differently
	'''
	z_dim, y_dim, x_dim = dimensions_3d
	return array_1d.reshape(x_dim,y_dim,z_dim)

def get_q_pixel_array(vti,name=None,threshold=True):
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


def get_umag_pixel_array(vti,name=None):
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



