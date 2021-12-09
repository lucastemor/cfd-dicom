"""
series writing, metadata management, etc 
"""

import numpy as np
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
import pydicom
import time,datetime,os,random


def check_and_convert_array(pixel_array):
	'''
	DICOM allows for custom bit depth UP TO 16 bit ...
	Should try and change this to something lower as we don't need all of that dynamic range
	13 bit, for example should be sufficient
	Note that this must also change alongside 'BitsAllocated', 'BitsStored','HighBit' tags
	'''
	if pixel_array.dtype != np.uint16:
		pixel_array = pixel_array.astype(np.uint16)
	return pixel_array.tobytes()

class dicom_stack:
	"""
	Based on metadata of working CT stack from Nicole
	Arbitrary coordinte system and measurements are used right now
	Can preserve physical measurements by instead encoding voxel length spacing for slice thickness, pixel spacing, and other tags related to spatial positions
	"""
	def __init__(self,write_dir,n_timesteps,case_name='case_name',patient_id='1234'):

		self.write_dir = write_dir
		self.n_timesteps = n_timesteps
		self.case_name  = case_name
		self.patient_id  = patient_id
		self.InStackPositionNumber = 0
		self.TemporalPositionIndex = 0
		self.LocalSliceId= 0

		self.date = datetime.datetime.now().strftime('%Y%m%d')
		self.studyID = str(int(random.random()*1e5))

		self.pixel_array = None
		self.slice_axis = None

		self.quantity = None

	def set_metadata(self):
		"""
		Initialize metadata
		Faster implementation is to inherit generic/shared metadata and only update the things that change in space or time
		"""

		header_metadata = {
		'MediaStorageSOPClass UID'      : '1.2.840.10008.5.1.4.1.1.2',
		'MediaStorageSOPInstanceUID'    : f'1.2.392.200036.9116.2.1220567263.{str(self.TemporalPositionIndex).zfill(10)}.{str(self.TemporalPositionIndex).zfill(6)}.{str(self.TemporalPositionIndex).zfill(3)}.{str(self.InStackPositionNumber).zfill(5)}',
		'ImplementationClassUID'        : "1.2.3.4",
		}

		filename = f'{self.write_dir}/IM-0001-{str(self.InStackPositionNumber).zfill(6)}.dcm'

		file_meta = FileMetaDataset()

		for entry in header_metadata:
			file_meta.__setattr__(entry,header_metadata[entry])

		ds = FileDataset(filename, {},
		                 file_meta=file_meta, preamble=b"\0" * 128)

		dimensions = self.pixel_array.shape

		#i.e., we're minimizing the number of files written at the expense of larger individual files (which can have internal compression)
		if self.slice_axis == 0:
			reconstruct_target = [self.LocalSliceId,0.0, 0.0]
			image_position_patient = pydicom.multival.MultiValue(pydicom.valuerep.DSfloat,[float(self.LocalSliceId),dimensions[0]/2.0,dimensions[1]/2.0])
		elif self.slice_axis ==1:
			reconstruct_target = [0.0,self.LocalSliceId, 0.0]
			image_position_patient = pydicom.multival.MultiValue(pydicom.valuerep.DSfloat,[dimensions[0]/2.0,float(self.LocalSliceId),dimensions[1]/2.0])
		elif self.slice_axis ==2:
			reconstruct_target = [0.0, 0.0, self.LocalSliceId]
			image_position_patient = pydicom.multival.MultiValue(pydicom.valuerep.DSfloat,[dimensions[0]/2.0,dimensions[1]/2.0,float(self.LocalSliceId)])

		file_metadata = {

		                 'ImageType' : ['DERIVED', 'SECONDARY', 'AXIAL', 'SUBTRACTION'],
		               'SOPClassUID' : '1.2.840.10008.5.1.4.1.1.2',
		            'SOPInstanceUID' : ds.file_meta.MediaStorageSOPInstanceUID,
		                 'StudyDate' : self.date,
		                'SeriesDate' : self.date,
		           'AcquisitionDate' : self.date,
		               'ContentDate' : self.date,		           
		                 'StudyTime' : '131420.000',
		                'SeriesTime' : '132018.033',
		           'AcquisitionTime' : str(132139.00 + 5*self.TemporalPositionIndex),
		               'ContentTime' : str(132139.00 + 5*self.TemporalPositionIndex),
		         'SeriesDescription' : f'CFD_{self.n_timesteps}_steps',
		         'SeriesInstanceUID' : f'1.3.12.2.1107.5.1.4.54693.3000000610250859320620000.{str(self.TemporalPositionIndex).zfill(4)}',
		                  'Modality' : 'OT',
		          'StudyDescription' : 'CFD',
		              'ProtocolName' : 'cfd-dicom',
		               'PatientName' : self.case_name,
		                 'PatientID' : self.patient_id,
		            'SliceThickness' : "1",
		               'TableHeight' : "0",
		         'RotationDirection' : 'CW',
		           'PatientPosition' : 'HFS',
		           'AcquisitionType' : 'STATIONARY',
 'ReconstructionTargetCenterPatient' : reconstruct_target,
		             'TablePosition' : self.LocalSliceId,
		          'StudyInstanceUID' : '1.3.6.1.4.1.12201.1011603639.75.20200110123223171.1',
		                   'StudyID' : self.studyID,
		              'SeriesNumber' : f"{self.n_timesteps}",
		            'InstanceNumber' : f"{self.InStackPositionNumber}",
		        'PatientOrientation' : ['L', 'P'],
		      'ImagePositionPatient' : image_position_patient,
		   'ImageOrientationPatient' : [1.00000, 0.00000, 0.00000, 0.00000, 1.00000, 0.00000],
		       'FrameOfReferenceUID' : '1.2.392.200036.9116.2.6.1.48.1220567263.1578629661.994251',
		             'SliceLocation' : f"{self.LocalSliceId}",
		                   'StackID' : f'1_{self.studyID}',
		    'In-StackPositionNumber' : self.InStackPositionNumber,
		     'TemporalPositionIndex' : self.TemporalPositionIndex,
		           'SamplesPerPixel' : 1,
		 'PhotometricInterpretation' : 'MONOCHROME2',
		                      'Rows' : dimensions[0],
		                   'Columns' : dimensions[1],
		              'PixelSpacing' : [1, 1],
		             'BitsAllocated' : 16,
	 	                'BitsStored' : 16,
		                   'HighBit' : 15,
		       'PixelRepresentation' : 1,
		              'WindowCenter' : "40.0",
		               'WindowWidth' : "80.0",
		          'RescaleIntercept' : "0.0",
		              'RescaleSlope' : "1.0",
		                 'PixelData' : check_and_convert_array(self.pixel_array),
		  'SmallestImagePixelValue'  : b'\\x00\\x00',
		    'LargestImagePixelValue' : b'\\xff\\xff',
		    	  'is_little_endian' : True,
		          'is_implicit_VR'   : True,

		}

		for entry in file_metadata:
			ds.__setattr__(entry,file_metadata[entry])

		#ds.SmallestImagePixelValue = b'\\x00\\x00'
		#ds.LargestImagePixelValue = b'\\xff\\xff'

		self.ds = ds


	def write_isotemporal_slices(self,array_3d):
		"""
		Given voxelized data at a given timestep, will write all slices, updating series_id
		Will advance timestep by one at the end
		Slice along smallest dimension 
		"""

		dimensions = array_3d.shape
		min_dimension = min(dimensions)
		slice_axis = dimensions.index(min_dimension)

		self.slice_axis = slice_axis

		for slice_id in range(min_dimension):
			slice_2d = np.take(array_3d,slice_id,slice_axis)
			self.pixel_array = slice_2d
			self.LocalSliceId = slice_id
			self.set_metadata()
			self.write()

			self.InStackPositionNumber += 1
		self.TemporalPositionIndex += 1


	def write(self):
		self.ds.save_as(self.ds.filename)


###############################
######### OLD #################
###############################

'''
class OLD_dicom_stack:
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
		"""
		For volumetric data
		"""
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

	def set_prerendered_image_metadata(self,pixel_array,timestep,angle):
		"""
		for frames rendered from paraview (or elsewhere)
		"""

		self.ds.Rows = pixel_array.shape[0]
		self.ds.Columns = pixel_array.shape[1]
		if pixel_array.dtype != np.uint16:
			pixel_array = pixel_array.astype(np.uint16)
		self.ds.PixelData = pixel_array.tobytes()

		self.ds.InstanceNumber = self.series_id 
		self.ds.ImagePositionPatient = [0,0,angle] 
		self.ds.SliceLocation = angle 

		self.ds.AcquisitionTime =  str(time.time())#str(timestep)
		self.ds.SeriesTime = str(float(timestep))
		self.ds.SeriesDescription = f'Timestep {timestep}'
		self.ds.SeriesNumber =  self.n_timesteps
		self.ds.SeriesInstanceUID = f'1.3.12.2.1107.5.1.4.54693.30000006102508593206200{str(int(timestep)).zfill(6)}'#f'1.3.12.2.1107.5.1.4.54693.3000000610250859320620000{str(int(timestep)).zfill(6)}'
		self.ds.FrameOfReferenceUID = f'1.3.12.2.1107.5.1.4.54693.30000006102508593206200{str(int(timestep)).zfill(6)}' #f'1.3.12.2.1107.5.1.4.54693.3000000610250859320620000{str(int(timestep)).zfill(6)}'

		self.ds.RotationDirection = 'CW'
		self.ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]


	def write_camera_positions(self,scaled_image,timestep,angle):
		"""
		for prerendered png frames, e.g., pathlines
		"""
		output_file = f'{self.write_dir}/IM-0001-{str(self.series_id).zfill(6)}.dcm'
		self.initialize_dicom(output_file)
		self.set_prerendered_image_metadata(scaled_image,timestep,angle)
		self.write()
		self.series_id+=1

	def write(self):
		self.ds.save_as(self.ds.filename)
'''
