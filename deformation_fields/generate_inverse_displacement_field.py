# Short script to generate INVERSE displacement field using Simple ITK in Python.
#
# v2025.08.18 - muniak@ohsu.edu

"""
First complete manual steps in FIJI (only done once, so didn't bother to automate...)

1)  Load each BrainReg deformation field TIFF into FIJI.
     Note: 'deformation_field_0', '_1', '_2' were renamed as 
           'ARA-to-CCFv3_25um_deformation_field_z', '_y', '_x' for clarity.

2)  Convert each TIFF into a _displacement_ field by:

    a) Scaling to pixel dimensions (from NiftyReg 1mm space to source voxel space, 
        e.g., *40 if 25um source voxels)... see: 
        https://forum.image.sc/t/whats-in-the-deformation-files/66409 and 
        https://github.com/brainglobe/cellfinder/blob/0dcba1c9fe125e37604735fb796cc9052d3ccf4b/cellfinder/analyse/analyse.py#L210-L226

    b) Subtracting z, y, or x coordinate from '_z', '_y', '_x', respectively.
    
    Thus, for each TIFF in FIJI, apply Process > Math > Macro... and use:
         "(v*40) - z" (or y or x)
    
    NOTE: The provided deformation field TIFFs have already had part (a) (*40) applied!!!

3)  Combine '_z', '_y', '_x' TIFFs as a single stack using 
     Image > Stacks > Tools > Concatenate... (uncheck "4D").
     
4)  Use Image > Properties... to set 25 um scale in all dimensions.

5)  Use File > Save As > MHD/MHA... to save as MHD/RAW for import to this python script.

NOTE: Further FIJI instructions at end of script!
"""

import os
import numpy as np
import SimpleITK as sitk

BASE_PATH = r'</PATH/TO/.MHD/FILE>'

# Read in RAW file output from FIJI.
original_dpfld = sitk.ReadImage(os.path.join(BASE_PATH, r'ARA-to-CCFv3_25um_displacement_field.mhd'), sitk.sitkFloat32)

# Raw file bytes are not ordered in the way sitk wants for building a DisplacementField.
# I can't figure out how to do so from FIJI (output is xyz"c"), so reordering here.
# There's probably a more efficient way to do this, but only had to be done once!
# Basically, convert to numpy array, reorder as 4D stack, then rebuild as sitkVectorField64 DisplacementField.
tmp = sitk.GetArrayFromImage(original_dpfld)
z_size = tmp.shape[0] // 3
tmp = np.stack([tmp[z:(z+z_size)] for z in np.arange(0, tmp.shape[0], z_size)[::-1]], axis=3)  # Note reversal of stack order.
dpfld = sitk.GetImageFromArray(tmp.astype(np.float64), isVector=True)

# Create transform then inverse.
dt = sitk.DisplacementFieldTransform(dpfld)
dt.SetInterpolator(sitk.sitkLinear)
tmp = sitk.InvertDisplacementField(dt.GetDisplacementField())
dti = sitk.DisplacementFieldTransform(tmp)
dti.SetInterpolator(sitk.sitkLinear)
inv_dpfld = dti.GetDisplacementField()

# Get inverse as array to resplit "channels", and save each as TIF stack.
tmp = sitk.GetArrayFromImage(inv_dpfld).astype(np.float32)

# Note re-reversal of stack order.
sitk.WriteImage(sitk.GetImageFromArray(tmp[:, :, :, 0].astype(np.float32), isVector=False), os.path.join(BASE_PATH, r'CCFv3-to-ARA_25um_displacement_field_x.tif'))
sitk.WriteImage(sitk.GetImageFromArray(tmp[:, :, :, 1].astype(np.float32), isVector=False), os.path.join(BASE_PATH, r'CCFv3-to-ARA_25um_displacement_field_y.tif'))
sitk.WriteImage(sitk.GetImageFromArray(tmp[:, :, :, 2].astype(np.float32), isVector=False), os.path.join(BASE_PATH, r'CCFv3-to-ARA_25um_displacement_field_z.tif'))

"""
Final FIJI steps.

6)  Load displacement field TIFFs into FIJI.
     'CCFv3-to-ARA_25um_displacement_field_z', '_y', '_x'

7)  Convert each TIFF back into a _deformation_ field by adding z, y, or x coordinate.
    
    Thus, for each TIFF in FIJI, apply Process > Math > Macro... and use:
         "v + z" (or y or x)
     
8)  Use Image > Properties... to set 25 um scale in all dimensions.

9)  Save each TIFF as 'CCFv3-to-ARA_25um_deformation_field_z', '_y', '_x'

10) Delete intermediary files if desired, only need deformation fields.
"""