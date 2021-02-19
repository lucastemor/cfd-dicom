"""
commandline tool for automating conversions and i/o
"""

import vtk,sys
import numpy as np
import glob,re,os

import preprocessing
import write_dicom as wd

from skimage import io 
from skimage.color import rgb2gray 


def glob_sort_files(series_path,extension):
	return sorted(glob.glob(series_path+f'/*.{extension}'),key = lambda x : int(re.findall(r'([\d]+).{}'.format(extension),x)[0]))#[::5]

if __name__=='__main__':
	"""
	Series id tracks "global" posistion in the dicom stack
	Slice id tracks "local" position in the timestep stack
	"""
	vti_series_path = '/home/lucas/Documents/viz/renders/Horos/surgical/MCA07/q_200x0200'
	png_series_path = '/Users/lucas/Documents/School/BSL/cfd-dicom/foo/pathline/'
	case_name = 'mca_07_q'
	quantity = 'q' #'u' # or 'q' or 'pathline' or 'dye'

	outdir = f'./output/{case_name}/'
	if os.path.exists(outdir) == False:
		os.mkdir(outdir)


	name = None

	if quantity == 'u':
		conversion_function = preprocessing.get_umag_pixel_array
	elif quantity == 'q':
		conversion_function = preprocessing.get_q_pixel_array
	elif quantity == 'dye':
		conversion_function = preprocessing.get_q_pixel_array
		name = 'n_points_in_cell'
	elif quantity == 'pathline':
		conversion_function = None
	else:
		print (f'quantity {quantity} unknown or unspecified')
		sys.exit()

	if quantity != 'pathline':
		time_series_files = glob_sort_files(vti_series_path,'vti')
		dicom_stack = wd.dicom_stack(write_dir=outdir,n_timesteps=len(time_series_files),case_name=case_name)

		for t,timestep in enumerate(time_series_files):
			reader = vtk.vtkXMLImageDataReader()
			reader.SetFileName(timestep)
			reader.Update()
			vti = reader.GetOutput()
		
			array_3d = conversion_function(vti,name=name)
			print (f'-------------Step {t}-------------------')
			print (array_3d.min(),array_3d.max())
			print (array_3d.astype(np.uint16).min(),array_3d.astype(np.uint16).max())
			dicom_stack.write_isotemporal_slices(array_3d)

	else:
		time_series_files = glob_sort_files(png_series_path,'png')

		#make sure files are written properly such that these data can be properly accessed from the filename
		#must be set in paraview
		angles = sorted(list(set([float(re.findall(r'angle([\d]+)',image)[0]) for image in time_series_files])))
		times = sorted(list(set([float(re.findall(r'time([\d]+)',image)[0]) for image in time_series_files])))

		#this is a bit of a bottleneck -- maybe a faster way to do this or a way to avoid it completely
		global_max_pixel_value = max([rgb2gray(io.imread(i)).max() for i in time_series_files])

		dicom_stack = wd.dicom_stack(write_dir=outdir,n_timesteps=len(time_series_files))

		for t in range(len(times)):
			for a in angles:
				image = rgb2gray(io.imread(time_series_files[dicom_stack.series_id]))
				scaled_image = preprocessing.global_png_scale(image,global_max_pixel_value)
				dicom_stack.write_camera_positions(scaled_image,t,a) 



