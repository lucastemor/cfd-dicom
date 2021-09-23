'''
Interactive selection of ROI for voxelization -- return a bounding box

Adapted code from Dan!
https://github.com/demacdo/geometry_tools/blob/master/geometry_tools/common.py

'''

import numpy as np
import pyvista as pv

class BoxSelect():
	'''Interactively define clipping box for mesh'''
	


class SacSelectTool():
    """ Interactively mark points using a probe.
    """ 
    def __init__(self, surf):
        self.surf = surf
        self.surf.point_data['Mask'] = np.zeros(self.surf.n_points)
        self.surf.point_data['TempMask'] = np.zeros(self.surf.n_points)

    def mask(self, center, radius):
        sphere = pv.Sphere(radius=radius, center=center)
        self.surf = self.surf.select_enclosed_points(sphere)
        mask_index = self.surf.point_data['SelectedPoints']
        ids = [x for x in range(self.surf.n_points) if mask_index[x] == True]
        self.selection = self.surf.extract_points(ids)

        self.selection.point_data['vtkOGIds'] = self.selection.point_data['vtkOriginalPointIds'].copy()
        self.selection = self.selection.extract_largest()
        mask = self.selection.point_data['vtkOGIds']
        self.selection = pv.PolyData(self.selection.points, self.selection.cells)
        self.selection.point_data['vtkOGIds'] = mask
        self.selection = self.selection.clean()

        self.surf.point_data['TempMask'] = self.surf.point_data['Mask'].copy()
        self.surf.point_data['TempMask'][self.selection.point_data['vtkOGIds']] = 1
            
    def select(self):

        def sphere_cb(xyz, probe):
            # Select enclosed points
            self.sphere = pv.Sphere(radius=probe.GetRadius(), center=xyz)
            # self.probe = probe
            self.p.add_mesh(self.sphere, color='r', opacity=0.4, name='probe')
            self.center = probe.GetCenter()
            self.radius = probe.GetRadius()
            self.mask(self.center, self.radius)
            self.p.add_mesh(self.surf, scalars='TempMask', name='surf')

        def choose_cb():
            self.surf.point_data['Mask'] = self.surf.point_data['TempMask']
            
        self.p = pv.Plotter()
        self.p.add_mesh(self.surf, scalars='TempMask', opacity=1.0, name='surf')
        self.p.add_sphere_widget(
            callback=sphere_cb, 
            center=np.mean(self.surf.points, axis=0), 
            radius=4.0, 
            pass_widget=True,
            theta_resolution=8,
            phi_resolution=8,
            style='wireframe',
            )
        self.p.add_key_event('space', choose_cb)
        self.p.show()

        if np.all(self.surf.point_data['Mask'] == 0):
            self.surf.point_data['Mask'] = self.surf.point_data['TempMask']