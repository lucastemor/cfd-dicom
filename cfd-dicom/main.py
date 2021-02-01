"""
commandline tool for automating conversions and i/o
"""

import vtk,sys
import numpy as np
import glob,re,os

import vti_to_dicom as vd
import write_dicom as wd


def glob_sort_vti(vti_series_path):
	return sorted(glob.glob(vti_series_path+'/*.vti'),key = lambda x : int(re.findall(r'([\d]+).vti',x)[0]))#[::5]

if __name__=='__main__':
	"""
	Series id tracks "global" posistion in the dicom stack
	Slice id tracks "local" position in the timestep stack
	"""
	vti_series_path = '/Users/lucas/Documents/School/BSL/cfd-dicom/foo/'
	case_name = 'testing'
	quantity = 'u' # or 'q'

	outdir = f'./output/{case_name}/'
	if os.path.exists(outdir) == False:
		os.mkdir(outdir)

	if quantity == 'u':
		conversion_function = vd.get_umag_pixel_array
	elif quantity == 'q':
		conversion_function = vd.get_q_pixel_array
	else:
		print (f'quantity {quantity} unknown or unspecified')
		sys.exit()


	time_series_files = glob_sort_vti(vti_series_path)
	dicom_stack = wd.dicom_stack(write_dir=outdir,n_timesteps=len(time_series_files))

	for t,timestep in enumerate(time_series_files):
		reader = vtk.vtkXMLImageDataReader()
		reader.SetFileName(timestep)
		reader.Update()
		vti = reader.GetOutput()
	
		array_3d = conversion_function(vti)
		dicom_stack.write_isotemporal_slices(array_3d)

