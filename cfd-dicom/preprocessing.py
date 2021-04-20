import vtk,sys
import numpy as np
import glob,re,os

from scipy.signal import convolve2d as conv2

MAX_DICOM_VALUE = 2^16 - 1


def global_png_scale(image,global_max_pixel_value):
	return (image*(MAX_DICOM_VALUE/global_max_pixel_value) + MAX_DICOM_VALUE).astype(np.uint16)




def get_umag_pixel_array(vti,name=None):
	"""
	for convenience
	"""
	return rescale_to_dicom(get_umag_voxel_array(vti,name))

'''
def rescale_to_dicom(array):
	"""
	data need to be rescaled so the pixel values are visible in dicom viewer
	2047 is specified by default based on  example files but hasn't relaly been played with
	"""
	return array * (MAX_DICOM_VALUE/array.max()) + array
'''
def rescale_to_dicom(array):
	return array*1000

def get_q_pixel_array(vti,name=None,threshold=True):
	"""
	Default name should be Q-criterion from paraview, butcan specify custom name if needed
	Can be used for other scalar quantities if cusotm name is passed
	"""

	dimensions = vti.GetDimensions()

	if name==None:
		name = 'Q-criterion'

	q_criterion = np.array(vti.GetPointData().GetArray(name))
	if threshold: 
		q_criterion = np.where(q_criterion < 0, 0.02, q_criterion)
		print ('Threshold!', q_criterion.min(), q_criterion.max()) 	

	return rescale_to_dicom(q_criterion.reshape(dimensions))

def get_umag_voxel_array(vti,name=None):
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
	
	return u_mag.reshape(dimensions)


def psf_volume(array_3d):
	for slice_id,slice_2d in enumerate(array_3d):
		array_3d[slice_id] = psf_slice(slice_2d)
	return(array_3d)


def psf_slice(slice_2d):
	psf = np.ones((5, 5)) / 25
	blurred = conv2(slice_2d, psf, 'same')
	return(blurred)
	



