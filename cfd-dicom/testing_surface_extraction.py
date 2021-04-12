import pyvista as pv
import numpy as np
from scipy.interpolate import griddata
from scipy.spatial import distance 
import matplotlib.pyplot as plt

from pyevtk.hl import imageToVTK
from scipy.ndimage import gaussian_filter

def closest_node(point,nodes):
	"""
	see https://codereview.stackexchange.com/questions/28207/finding-the-closest-point-to-a-list-of-points
	"""
	closest_index = distance.cdist([point], nodes).argmin()
	return closest_index


wss_path = '/mnt/3414B51914B4DED4/dicom/data/surface_data/wss_surface.vtu'

#load surface
surface = pv.read(wss_path)

#slice along z-axis
n_slices = 200
slice_dim_x = 200
slice_dim_y = 200

voxel_array = np.zeros((slice_dim_x,slice_dim_y,n_slices))

#get bounding box
xmin, xmax, ymin, ymax, zmin, zmax = surface.bounds


x = np.linspace(xmin,xmax,num=slice_dim_x)
y = np.linspace(ymin,ymax,num=slice_dim_y)
z = np.linspace(zmin,zmax,num=n_slices)

grid_points = np.array([[xx,yy] for xx in x for yy in y])\

slices = surface.slice_along_axis(n=n_slices,axis='z')

for slice_id in range(n_slices):
	resampled_image = np.zeros(len(grid_points)) # will reshape into 2D image after computation
	slice_name = f'slice{slice_id}'
	model_slice = slices[slice_name]

	planar_points = model_slice.points[:,0:2]
	values = model_slice.point_arrays['wss_mag']

	for point_id,point in enumerate(planar_points):
		closest_node_id = closest_node(point,grid_points)
		resampled_image[closest_node_id] = values[point_id]

	voxel_array[:,:,slice_id] = resampled_image.reshape(slice_dim_x,slice_dim_y)




#see https://vtk.org/Wiki/VTK/Writing_VTK_files_using_python
#dimensions
nx, ny, nz = voxel_array.shape
ncells = nx * ny * nz
npoints = (nx + 1) * (ny + 1) * (nz + 1) 
	
outfile = f'/mnt/3414B51914B4DED4/dicom/data/surface_data/200x200x200sigma2wss_subtraction'

spacing = tuple(np.abs([(xmax-xmin)/slice_dim_x,(ymax-ymin)/slice_dim_y,(zmax-zmin)/n_slices]))


blurred = 3*gaussian_filter(voxel_array,sigma=2)

voxel_array += np.where(voxel_array==0,blurred,0)

imageToVTK(outfile, origin= (xmin,ymin,zmin), spacing=spacing, pointData = {"wss_mag" : voxel_array} ) #can work for cell data too


#resampled = griddata(planar_points,values,grid_points,'nearest').reshape(slice_dim_x,slice_dim_y)
#plt.pcolormesh(resampled_image)
#plt.show()

#slices['slice0'].point_arrays['wss_mag']

#resample slice onto grid
#resamp = grid.sample(slice0)

#vox = pv.voxelize(surface,check_surface=False, density=0.1)


#p = pv.Plotter()
#p.add_mesh(vox)
#p.show()