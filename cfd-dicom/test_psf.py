import vtk,sys
import numpy as np
import glob,re,os

import preprocessing
import write_dicom as wd

from skimage import io 
from skimage.color import rgb2gray 

import matplotlib.pyplot as plt

from scipy.signal import convolve2d as conv2
	
vti_path = '/Users/lucas/Documents/School/BSL/Horos/voxelized_pathline/grid200cube_mask200_track25_radius0p05_pathline_150.vti'


reader = vtk.vtkXMLImageDataReader()
reader.SetFileName(vti_path)
reader.Update()
vti = reader.GetOutput()


array_3d = preprocessing.get_umag_pixel_array(vti)

slice_id = 100

image = array_3d[slice_id]

plt.imshow(image, cmap='gray')   

psf = np.ones((5, 5)) / 25
blurred = conv2(image, psf, 'same')


#fig,ax = plt.subplots(1,2)



plt.imshow(blurred, cmap='gray') 