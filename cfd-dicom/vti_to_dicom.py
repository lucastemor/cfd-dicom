import vtk,sys
import numpy as np
import glob,re,os

import write_dicoms

def glob_sort_vti(vti_series_path):
	return sorted(glob.glob(vti_series_path+'/*.vti'),key = lambda x : int(re.findall(r'([\d]+).vti',x)[0]))#[::5]

def get_q_pixel_array(vti,name=None,threshold=False,max_dicom_value=2047):
	"""
	for convenience
	"""
	return rescale_to_dicom(get_q_voxel_array(vti,name,threshold),max_dicom_value)

def get_umag_pixel_array(vti,name=None, max_dicom_value=2047):
	"""
	for convenience
	"""
	return rescale_to_dicom(get_umag_voxel_array(vti,name),max_dicom_value)

def rescale_to_dicom(array, max_dicom_value=2047):
	"""
	data need to be rescaled so the pixel values are visible in dicom viewer
	2047 is specified by default based on  example files but hasn't relaly been played with
	"""
	return array * (max_dicom_value/array.max()) + array

def get_q_voxel_array(vti,name=None,threshold=False):
	"""
	Default name should be Q-criterion from paraview, butcan specify custom name if needed
	Can be used for other scalar quantities if cusotm name is passed
	"""

	dimensions = vti.GetDimensions()

	if name==None:
		name = 'Q-criterion'

	q_criterion = np.array(vti.GetPointData().GetArray('Q-criterion'))
	if threshold: q_criterion = np.where(q_criterion < 0, 0, q_criterion) #double check thresholding 	

	return q_criterion.reshape(dimensions)

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

#def array_3d_to_dicom_series()


if __name__=='__main__':
	vti_series_path = '/Users/lucas/Documents/School/BSL/cfd-dicom/foo/'
	case_name = 'testing'
	quantity = 'u' # or 'q'

	#if case_name not in os.
		#os.mkdir

	if quantity == 'u':
		conversion_function = get_umag_pixel_array
	elif quantity == 'q':
		conversion_function = get_q_pixel_array
	else:
		print (f'quantity {quantity} unknown or unspecified')
		sys.exit()


	time_series_files = glob_sort_vti(vti_series_path)
	series_id = 0
	for t,timestep in enumerate(time_series_files):
		reader = vtk.vtkXMLImageDataReader()
		reader.SetFileName(timestep)
		reader.Update()
		vti = reader.GetOutput()
	
		array_3d = conversion_function(vti)

		for slice_id, slice_array_2d in enumerate(array_3d):
			plane_centre_coords = pydicom.multival.MultiValue(pydicom.valuerep.DSfloat,[dimensions[0]/2.0,dimensions[1]/2.0,float(i)]) #each slice 2d is an xy plane slice thgat moves along z axis





