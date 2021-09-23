# cfd-dicom
Tools for converting CFD data to 4D dicom series

Args are specified at the top of the `if __name__ == '__main__'` block.

##important: need to fix coord. system transformation for pathlines

### Usage - `cfd-dicom.py` args
**Scalars**
1. **`mesh_path`** - path to `.h5` mesh containing `coordinates` and `topology`
2. **`n_steps`*** - number of timesteps to sample. Will pick **uniformly** over all avaliable timesteps. E.g., if there are 1000 steps total and `n_steps = 100` every 10th step will be written to DICOM
3. **`voxel_size`** - isotropic voxel length, assumed to be in the same length scale as the unstructrued mesh, usually mm 
4. **`quantities`** - velocity, q-criterion, or pathlines. Can specify any combination of `u`, `q`, `path` depending on the DICOM series to be generated. **Note**: if `path` is specified, pathline args shown beow must also be specified

**Pathlines**
5. **`particle_polydata_path`** - path to `.vtp` series containing tracked particle polydata in the same mesh as specified above. Assumes data are tracked particles written out of Paraview's `ParticleTracer` filter 
6. **`sigma`** - Radius of blurring amount. Assumes same length scales as mesh/voxel size above.
7. **`track_length`** - Akin to "exposure time" (Steinman, 2000). Larger values will result in longer particle streaks. A `track_length` of `n` means that the particle field is imaged every `n` timesteps. **Note**: In this code there is a difference between "imaging" the flow (i.e., resetting the "light" on the lens) and actually *writing* a DICOM snapshot. The former is dependent on `track_length` and the later on  `n_steps`. What is ultimately written to DICOM is the *intersection* of frames accounted for in the sampling with `n_steps` and the snapshot frames resulting from `track_length`
8. **`step_stride`** - Steps over which to compute the pathlines. Akin to frame rate from Steinman, 2000. For example, if `step_stride=2`, pathlines will be "drawn" between tracked particles at `t=0, t=2, t=4, etc` as opposed to `step_stride=1` where lines will be drawn between each successive timestep (i.e., `t=0, t=1, t=2, t=3, etc`)

Once these are specified you may run `cfd-dicom.py` and will have to specify a bounding box for the ROI you wish to voxelize.

After running `cfd-dicom.py` DICOM series will be in `./output`

The most common errors will probably result from inconsistent naming conventions of the `.h5` mesh or time series files and tracked `.vtp` particles. 

Some process work and/or other interesting (but no longer necessary) scripts can be founf in `./extra`

### Dependencies
`pyvista`, `h5py`, `numpy`, `scipy`, `pathlib`, `skimage`, `pydicom`
