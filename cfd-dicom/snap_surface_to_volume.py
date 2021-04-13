import pyvista as pv
import numpy as np
import vtk
from scipy.interpolate import griddata
from scipy.spatial import distance 
import matplotlib.pyplot as plt

from pyevtk.hl import imageToVTK
from scipy.ndimage import gaussian_filter

import write_dicom as wd

def point_id_to_structured_coords(vti, point_id):
	"""
	vti written in paraview will count from (0,0,0) to (dim,dim,dim) incrementing xdime+1 every time and periodically loopiung back to x=0 but then iuncrmeenting y up to y dim .. etc..
	"""

	x_dim,y_dim,z_dim = vti.GetDimensions()

	x = point_id%x_dim
	y = int(((point_id-x)/x_dim)%y_dim)
	z = int((point_id-x-y*y_dim)/(x_dim*y_dim))

	return x,y,z



wss_path = '/mnt/3414B51914B4DED4/dicom/data/surface_data/wss_surface.vtu'
vti_path = '/home/lucas/Documents/viz/renders/Horos/surgical/MCA07/velocity_200x200/velocity200x200surgical_mca07_0.vti'

#load surface
surface = pv.read(wss_path)

#load voxelized volume

reader = vtk.vtkXMLImageDataReader()
reader.SetFileName(vti_path)
reader.Update()
vti = reader.GetOutput()


voxel_array = np.zeros(vti.GetDimensions())
point_count = np.zeros_like(voxel_array)

#instea of slicing as before, we move each point to closest unsturctured point, and where there are overlaps, we average
for point,wss in zip(surface.points,surface.point_arrays['wss_mag']):
	point_id_on_grid = vti.FindPoint(point)

	if point_id_on_grid != -1:
		x,y,z = point_id_to_structured_coords(vti, point_id_on_grid)
		voxel_array[x,y,z] += wss
		point_count[x,y,z] += 1

#add 999s to remaining zero indices so averaging doesnt have 0 divisor - won't effect anything becasue WSS should be zero at these points anyways
point_count = np.where(point_count==0,999,point_count)

#average wss where there is overlap
voxel_array = voxel_array/point_count




outdir = '/mnt/3414B51914B4DED4/dicom/data/surface_data/DICOM_1'
case_name = 'test_wss_snappoints'
n_timesteps= 1
dicom_stack = wd.dicom_stack(write_dir=outdir,n_timesteps=n_timesteps,case_name=case_name)

voxel_array *= 1000
dicom_stack.write_isotemporal_slices(voxel_array)

'''
#see https://vtk.org/Wiki/VTK/Writing_VTK_files_using_python
#dimensions
nx, ny, nz = voxel_array.shape
ncells = nx * ny * nz
npoints = (nx + 1) * (ny + 1) * (nz + 1) 
	
outfile = f'/mnt/3414B51914B4DED4/dicom/data/surface_data/200x200x200_snaptogrid_blurred'


#blurred = 3*gaussian_filter(voxel_array,sigma=2)

#voxel_array += np.where(voxel_array==0,blurred,0)

imageToVTK(outfile, origin= vti.GetOrigin(), spacing=vti.GetSpacing(), pointData = {"wss_mag" : voxel_array} ) #can work for cell data too
'''

#resampled = griddata(planar_points,values,grid_points,'nearest').reshape(slice_dim_x,slice_dim_y)
#plt.pcolormesh(resampled_image)
#plt.show()

#slices['slice0'].point_arrays['wss_mag']

#resample slice onto grid
#resamp = grid.sample(slice0)
