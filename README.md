# cfd-dicom
Tools for converting CFD data to 4D dicom series

At its most basic level, you can use `cfd-dicom.py` to convert scalars or integrated quantities to DICOM series to be rendered in Horos. Input data must be specified in the `cfd-dicom.py` script. In the future, could be modified to do batch jobs or pass args as command line args. Ultimately, DICOM series are stored to `./cfd-dicom/output/{case folder}_voxelsize={voxel_size}_nsteps={n_steps}/`. Descriptions of `voxel_size`, `n_steps`, and all other necesary inputs are given below. 

### Usage - `cfd-dicom.py` args. 

For now, these are immeadiately following the `if __name__ == '__main__'` line.

By default, these are specified to read the example data in `cfd-dicom/example_data/`

**Scalars**
1. **`mesh_path`** - path to `.h5` mesh containing `coordinates` and `topology`
2. **`n_steps`*** - number of timesteps to sample. Will pick **uniformly** over all avaliable timesteps. E.g., if there are 1000 steps total and `n_steps = 100` every 10th step will be written to DICOM. Eventually, more intelligent sampling schemes could be used, like "bulllet time"
3. **`voxel_size`** - isotropic voxel length, assumed to be in the same length scale as the unstructrued mesh, usually mm 
4. **`quantities`** - velocity, q-criterion, or pathlines. Must be specified as a list. Can specify any combination of `u`, `q`, `path` depending on the DICOM series to be generated. **Note**: if `path` is specified, pathline args shown below must also be specified. If you are not interested in pathlines or streamlines, you can stop here and run the script.

**Pathlines**

5. **`particle_polydata_path`** - path to `.vtp` series containing tracked particle polydata in the same mesh as specified above. Assumes data are tracked particles written out of Paraview's `ParticleTracer` filter. If you use a naming convention different from Praview's defaults, e.g., `ParticleId` for particle ids, you will have to change indexing in relevant parts of the code.  
6. **`sigma`** - Radius of blurring amount. Assumes same length scales as mesh/voxel size above.
7. **`track_length`** - Akin to "exposure time" (Steinman, 2000). Larger values will result in longer particle streaks. A `track_length` of `n` means that the particle field is imaged every `n` timesteps. **Note**: In this code there is a difference between "imaging" the flow (i.e., resetting the "light" on the lens) and actually *writing* a DICOM snapshot. The former is dependent on `track_length` and the later on  `n_steps`. What is ultimately written to DICOM is the *intersection* of frames accounted for in the sampling with `n_steps` and the snapshot frames resulting from `track_length`
8. **`step_stride`** - Steps over which to compute the pathlines. Akin to frame rate from Steinman, 2000. For example, if `step_stride=2`, pathlines will be "drawn" between tracked particles at `t=0, t=2, t=4, etc` as opposed to `step_stride=1` where lines will be drawn between each successive timestep (i.e., `t=0, t=1, t=2, t=3, etc`)

Once these are specified you may run `cfd-dicom.py` and will have to specify a bounding box for the ROI you wish to voxelize.

The most common errors will probably result from inconsistent naming conventions of the `.h5` mesh or time series files and tracked `.vtp` particles. E.g., sometimes veocity is saved as `u` (newer) sometimes as `Solution/u` (older)

Some process work and/or other interesting (but no longer necessary) scripts can be founf in `./extra`

### Dependencies
`pyvista`, `h5py`, `numpy`, `scipy`, `pathlib`, `skimage`, `pydicom`


### To do
1. Streamlines and/or velocity interpolation
2. There are some redundant/repetitive chunks of code, there is probably a better, more "streamlined ;-)" approach to organization 