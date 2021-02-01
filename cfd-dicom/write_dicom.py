"""
series writing, metadata management, etc 
"""

import numpy as np
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
import pydicom
import time,datetime,os


class dicom_stack:
	"""
	Main class for setting metadata and writing in 4D
	Kind of like a wrapper for pydicom dataset object 

	A lot of this is pulled form this link below and hasn't really been explored completely

	Useful resource: https://stackoverflow.com/questions/14350675/create-pydicom-file-from-numpy-array
	Shared/common metadata for either volumetric data or prerendred images

	"""
	def __init__(self,write_dir,n_timesteps,case_name='case_name',patient_id='1234'):

		self.write_dir = write_dir
		self.n_timesteps = n_timesteps
		self.case_name  = case_name
		self.patient_id  = patient_id
		self.series_id = 0
		self.timestep = 0

	def write_isotemporal_slices(self,array_3d):
		"""
		Given voxelized data at a given timestep, will write all slices, updating series_id
		Will advance timestep by one at the end
		"""
		dimensions = array_3d.shape

		for slice_id,slice_2d in enumerate(array_3d):
			output_file = f'{self.write_dir}/IM-0001-{str(self.series_id).zfill(6)}.dcm'
			self.plane_centre_coords = pydicom.multival.MultiValue(pydicom.valuerep.DSfloat,[dimensions[0]/2.0,dimensions[1]/2.0,float(slice_id)]) #each slice 2d is an xy plane slice thgat moves along z axis	
			self.slice_id = slice_id
			self.initialize_dicom(output_file)
			self.set_volumetric_slice_metadata(slice_2d)
			self.write()

			self.series_id += 1
		self.timestep += 1


	def initialize_dicom(self,filename):
		"""
		dicom header etc.. A lot of this is pulled form this link below and hasn't really been explored completely

		Useful resource: https://stackoverflow.com/questions/14350675/create-pydicom-file-from-numpy-array
		Shared/common metadata for either volumetric data or prerendred images
		"""

		filename_little_endian = filename

		file_meta = FileMetaDataset()
		file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.2'
		file_meta.MediaStorageSOPInstanceUID = "1.2.3"
		file_meta.ImplementationClassUID = "1.2.3.4"

		ds = FileDataset(filename_little_endian, {},
		                 file_meta=file_meta, preamble=b"\0" * 128)

		ds.PatientName = self.case_name
		ds.PatientID = self.patient_id

		# Set the transfer syntax
		ds.is_little_endian = True
		ds.is_implicit_VR = True

		# Set creation date/time
		dt = datetime.datetime.now()
		ds.ContentDate = dt.strftime('%Y%m%d')
		timeStr = dt.strftime('%H%M%S.%f')  # long format with micro seconds
		ds.ContentTime = timeStr

		ds.SeriesNumber = self.n_timesteps

		ds.SamplesPerPixel = 1
		ds.PhotometricInterpretation = "MONOCHROME2"
		ds.PixelRepresentation = 0
		ds.HighBit = 15
		ds.BitsStored = 16
		ds.BitsAllocated = 16
		ds.SmallestImagePixelValue = b'\\x00\\x00'
		ds.LargestImagePixelValue = b'\\xff\\xff'

		self.ds = ds

	def set_volumetric_slice_metadata(self,pixel_array):
		self.ds.Columns = pixel_array.shape[0]
		self.ds.Rows = pixel_array.shape[1]
		if pixel_array.dtype != np.uint16:
			pixel_array = pixel_array.astype(np.uint16)
		self.ds.PixelData = pixel_array.tobytes()

		self.ds.InstanceNumber = self.slice_id
		self.ds.ImagePositionPatient = self.plane_centre_coords
		self.ds.SliceLocation = self.slice_id

		self.ds.AcquisitionTime =  str(time.time())#str(timestep)
		self.ds.SeriesTime = str(float(self.timestep))
		self.ds.SeriesDescription = f'Timestep {self.timestep}'
		self.ds.SeriesNumber =  self.n_timesteps #number of timesteps
		self.ds.SeriesInstanceUID = f'1.3.12.2.1107.5.1.4.54693.3000000610250859320620000{str(self.timestep).zfill(4)}'
		self.ds.file_meta.MediaStorageSOPInstanceUID = f"1.3.12.2.1107.5.1.4.54693.30000006101906583670300{str(self.series_id).zfill(6)}"

	def write(self):
		self.ds.save_as(self.ds.filename)



