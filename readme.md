# Striatum Clusters to CCFv3

**Functions to query striatum cluster voxels from Hunnicutt et al. (2016) -- Fig 5d.**   

Based on https://elifesciences.org/articles/19103  
Uses input coordinates in CCFv3 um space [AP, DV, LM].  
Returns cluster id values (1-15) based on top-to-bottom order of clusters within dendrogram shown in Fig 5c.  

v2025.10.13  
muniak@ohsu.edu  
Mao Lab - Vollum Institute

These tools were developed in python/napari, and then ported to Matlab.  Both versions of the code are self-contained (no other custom modules), but do require support files  contained in the sub-folders of this package.

Code related to the original 2016 paper can be located at:
https://github.com/BJHunnicutt/anatomy/  

napari_functions.py
-------------------
This script contains all code necessary to query the clusters.  The top section contains the core functions, and also sets up some visualizations in napari.  The bottom portion contains some example demonstration code.

query_striatum_map.m
--------------------
This is a self-contained Matlab function that allows you query the clusters in Matlab.  Use "help" in Matlab for info on usage.

atlas_templates/
----------------
25um coronal templates of CCFv3 and ARA atlases.  These images were used to calculate the mapping between ARA<->CCFv3 using brainreg.

csvs/
-----
Support files containing:
- Cluster ids/RGB colors for building colormaps.
- Scalings for data volumes at different resolutions to um space.
- Translations for striatum-specific data volumes at different resolutions to um space.

deformation_fields/
-------------------
Deformation fields for ARA>CCFv3 mapping, as calculated using brainreg, with a post-hoc scale adjustment (see .py file).  Also contains inverse deformation fields generated following the script *generate_inverse_displacement_field.py*.  The .json file  preserves the parameters used with brainreg.

striatum_clusters_ara/
----------------------
The original striatum cluster volumes as calculated in Hunnicutt et al. (2016).  Spacing is 150x150x150 um [AP, DV, LM].  Volumes can be aligned to ARA um space using the scale/translation CSVs, see napari scripting for details.

striatum_clusters_ccfv3/
------------------------
The result of mapping striatum clusters from ARA>CCFv3.  Original cluster voxels in ARA space were re-created as 3D "cubes", and the vertices of these cubes were mapped to CCFv3 based on the deformation fields.  Each voxel of CCFv3 at 25um resolution was tested to see if it was "inside" a mapped cube, and if so, that voxel was assigned the id value of the source ARA voxel from which the cube was derived.  Accordingly, these volumes exist at 25x25x25 um [AP, DV, LM] resolution in order to capture the 3D deformations due to mapping.  The code to generate these volumes is found at the bottom of *napari_functions.py* in the commented-out section.

striatum_masks/
------------------------
Two binary masks of the striatum:

*striatum_mask_ara-bjh_100um.tif* is based on the original 25x25x25 um mask annotated by Hunnicutt in ARA space (available at the original Github link above), which was subsequently downsampled to 100x100x100 um and from which a separate mask of the internal capsule was subtracted (only available at 100x100x100 um, see same Github code), giving this result.  This mask was later downsampled to 150x150x150 um to serve as the template for the cluster maps.  

*striatum_mask_ccf_25um.tif* is based on the CCFv3 atlas, combining the masks for dorsal and ventral striatum ('STRd' and 'STRv') and subtracting the mask for olfactory tubercle ('OT').  The code for this can be found at the bottom of *napari_functions.py* in the commented-out section.

napari env
----------
A working napari environment can be recreated by following the standard installation instructions:
https://napari.org/stable/tutorials/fundamentals/installation.html

To run generate_inverse_displacement_field.py, it is necessary to install SimpleITK, e.g., "conda install simpleitk".

To replicate the template mapping between ARA and CCFv3, install brainreg:
https://brainglobe.info/documentation/brainreg/installation.html#installation