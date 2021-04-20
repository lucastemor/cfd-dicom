import numpy as np
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
import pydicom, glob, difflib

def compare_metadata(image1,image2):
	'''
	https://pydicom.github.io/pydicom/stable/auto_examples/plot_dicom_difference.html
	'''

	datasets = (image1,image2)

	rep = []
	for dataset in datasets:
		lines = str(dataset).split("\n")
		lines = [line + "\n" for line in lines]  # add the newline to end
		rep.append(lines)


	diff = difflib.Differ()
	for line in diff.compare(rep[0], rep[1]):
		if line[0] != "?":
			print(line)




#orking_path_ts1 = '/home/lucas/Downloads/4DCTA (2 time steps)/Timestep 1/_SUB_Axial_SUB_Head_05_CE_5'
#working_path_ts2 = '/home/lucas/Downloads/4DCTA (2 time steps)/Timestep 2/_SUB_Axial_SUB_Head_05_CE_5'

#cfd_path_ts1 = '/home/lucas/Downloads/MCA07_vel_0_1000'

working_step1 = pydicom.read_file(sorted(glob.glob(working_path_ts1+'/*.dcm'))[1])
#working_step2 = pydicom.read_file(sorted(glob.glob(working_path_ts1+'/*.dcm'))[1])
working_step2 = pydicom.read_file(sorted(glob.glob(working_path_ts2+'/*.dcm'))[1])

cfd_step1 = pydicom.read_file(sorted(glob.glob(cfd_path_ts1+'/*.dcm'))[0])

image1 = pydicom.read_file('/home/lucas/Downloads/WORKING_ORIG.dcm')
image2 = pydicom.read_file('/home/lucas/Downloads/WORKING_EXP.dcm')
