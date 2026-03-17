# napari_functions.py
# v2025.10.13
# muniak@ohsu.edu
#
# Functions to query striatum cluster map(s) from Hunnicutt et al. (2016) Fig 5d. 
# Uses input coordinates in CCFv3 um space [AP, DV, LM].
# Returns cluster id values (1-15) based on top-to-bottom order of clusters within dendrogram shown in Fig 5c.

import os
import pandas as pd
from scipy.interpolate import RegularGridInterpolator as rgi
from scipy.spatial import distance
from skimage.io import imread

# Define local path containing this script and support files.
BASE_PATH =  r'</PATH/TO/DIRECTORY/CONTAINING/SCRIPT>'
# Cluster voxel variants.
CLUSTER_NUMS = [2, 3, 4, 15]

viewer = napari.current_viewer()

# Load scale/translation dicts.
vx_s = pd.read_csv(os.path.join(BASE_PATH, r'csvs/voxel_scales.csv'), index_col=0, header=0).to_dict('list')
vx_t = pd.read_csv(os.path.join(BASE_PATH, r'csvs/voxel_translations.csv'), index_col=0, header=0).to_dict('list')

# Load ARA template.
ara_im = imread(os.path.join(BASE_PATH, r'atlas_templates/ARA_25um_coronal_template.tif'))
ara = viewer.add_image(ara_im, name='ARA', scale=vx_s['25'], colormap='gray', visible=False)

# Load CCFv3 template.
ccf_im = imread(os.path.join(BASE_PATH, r'atlas_templates/CCFv3_25um_coronal_template.tif'))
ccf = viewer.add_image(ccf_im, name='CCFv3', scale=vx_s['25'], colormap='gray', visible=True)

# Load striatum masks.
mask_ara = imread(os.path.join(BASE_PATH, r'striatum_masks/striatum_mask_ara-bjh_100um.tif'))
viewer.add_image(mask_ara, name='striatum [ara-bjh]', scale=vx_s['100'], colormap='red', blending='additive', contrast_limits = [0, 1], visible=False)
mask_ccf = imread(os.path.join(BASE_PATH, r'striatum_masks/striatum_mask_ccf_25um.tif'))
viewer.add_image(mask_ccf, name='striatum [ccf]', scale=vx_s['25'], colormap='red', blending='additive', contrast_limits = [0, 1], visible=False)

# Load cluster color CSV.
df = pd.read_csv(os.path.join(BASE_PATH, r'csvs/cluster_colors.csv'), index_col=0, header=0)
# Convert to array and add row for empty space.
colors = np.vstack(( np.array([[0, 0, 0]]),
                     df.loc[:, ['r', 'g', 'b']], )).astype(float)
# Convert to 0.0 - 1.0 scale.
colors /= 255
# Add alpha...
colors = np.hstack((colors, np.ones((colors.shape[0], 1))))
# ...but set empty to transparent.
colors[0, -1] = 0
# Define colormap dict.
cluster_colormap = {
    'colors': colors,
    'name': 'striatum_clusters',
}

# Load original striatum cluster voxel TIFs.
clusters_ara = {c: imread(os.path.join(BASE_PATH, r'striatum_clusters_ara/striatum_{:02d}_clusters_ara.tif'.format(c))) for c in CLUSTER_NUMS}
for c in clusters_ara:
    viewer.add_image( clusters_ara[c], name='cluster[ara] {:02d}'.format(c),
                      scale=vx_s['150'], translate=vx_t['150'],
                      colormap=cluster_colormap, contrast_limits=[0, 15], visible=False )

# Load mapped striatum cluster voxel TIFs in 25um CCFv3 space.
clusters_ccf = {c: imread(os.path.join(BASE_PATH, r'striatum_clusters_ccfv3/striatum_{:02d}_clusters_ccf.tif'.format(c))) for c in CLUSTER_NUMS}
for c in clusters_ccf:
    viewer.add_image( clusters_ccf[c], name='cluster[ccf] {:02d}'.format(c),
                      scale=vx_s['25'], translate=None,
                      colormap=cluster_colormap, contrast_limits=[0, 15], visible=False )


def flip_over_midline(
        pts,  # N x 3 array of points.
    ):
    """ Flip any left side points over the midline to right side. (M-L == 5700 µm; dim == 2)
    """
    flp_pts = pts.copy()
    idx = pts[:, 2] < 5700
    flp_pts[idx, :] = ( flp_pts[idx, :] * np.array([1, 1, -1]) ) + np.array([0, 0, 11400.])
    return flp_pts
    

def show_voxels_as_points(
        c_vol,                           # Image containing cluster values.
        c_offset = np.array([0, 0, 0]),  # Offset of c_vol (in atlas um resolution).
        c_scale = np.array([1, 1, 1]),   # Scaling of c_vol.
        name = 'cluster as points',      # Custom name for point layer in napari.
        size = 5,                        # Sizing of points in napari.
        tform = None,                    # Set to 'ara' to map ARA voxels to CCFv3 space.
    ):
    """ Show voxel image as points in napari.
    """
    # Get indices of voxel centers.
    z, y, x = np.where(c_vol > 0)
    # Get values of voxels.
    v = c_vol[z, y, x]
    # Arrange indices as N x 3 array, scale to atlas resolution, and apply offset.
    zyx = ( np.array([z, y, x]).T * c_scale ) + c_offset
    # Map points, if requested.
    if tform:
        zyx = map_coordinates(zyx, df_id='ara')
    # Visualize points.
    # Dividing and scaling points image resolution to keep napari scrollbar nice.
    viewer.add_points(zyx / vx_s['25'], name=name, scale=vx_s['25'], size=5,
        face_color=cluster_colormap['colors'][v, :])


def get_cluster_assignments(
        pts,                             # Points to test in atlas um resolution [AP, DV, LM].
        c_vol,                           # Image containing cluster values.
        c_offset = np.array([0, 0, 0]),  # Offset of c_vol relative to pts, if any (in atlas um resolution).
        c_scale = np.array([1, 1, 1]),   # Scaling of c_vol relative to pts.
        map_to_nearest = False,          # For pts that are not within defined voxels, return value of nearest voxel.
                                         # This option also returns boolean array indicating which pts are 
                                         # _inside_ voxels (True) vs. outside/mapped (False).
                                         # Note: This option scales poorly with number of points for 'forward', sorry.
        return_distance = False,         # Return euclidean distances from pts to voxel centers.
        return_voxels = False,           # Return coordinates of matched voxels (in atlas um resolution).
    ):
    """ Returns cluster voxel assignments.
        If not mapping to nearest, unassigned/out-of-bounds voxels will have cluster value of 0.
    """
    # Init result array.
    v = np.zeros(pts.shape[0], dtype=int)
    # Adjust points for offset of cluster image.
    pts = pts - c_offset
    # If requested, compute distance between pts and voxel centers before rounding to indices.
    if return_distance or map_to_nearest:
        # Get indices of voxel centers.
        cz, cy, cx = np.where(c_vol > 0)
        # Arrange indices as N x 3 array and scale to atlas resolution.
        czyx = np.array([cz, cy, cx]).T * c_scale
        # Find euclidean distances between all points and voxel coordinates.
        pdists = distance.cdist(pts, czyx, 'euclidean')
        # Find nearest voxel coordinate for each test point.
        pdists_min_idx = pdists.argmin(axis=1)
        # Get euclidean distance for nearest voxel coordinate for each test point.
        pdists_min = pdists[range(len(pts)), pdists_min_idx]
    # First, find values for points that are contained in voxels, but not others.    
    # Scale points to sampling space of cluster image.
    pts = pts / c_scale
    # Resample as integers for indexing.
    pts = np.round(pts).astype(int)
    # Identify out-of-bounds points (automatically v==0).
    oob = np.any((pts >= np.array(c_vol.shape)) | (pts < 0), axis=1)
    # Return cluster values for in-bounds points only.
    v[~oob] = c_vol[pts[~oob, 0], pts[~oob, 1], pts[~oob, 2]]
    # Assemble atlas coordinates of assigned voxels.
    v_pts = ( pts.astype(float) * c_scale ) + c_offset
    # Get indices of unassigned points.
    v_null = v==0
    # Set unassigned voxels to nan.
    v_pts[v_null, :] = np.nan
    # If requested, return value of nearest voxel, even if outside of cluster map.
    if map_to_nearest:
        # Get cluster voxel coordinate that was nearest each test point.
        pts2 = czyx[pdists_min_idx, ...]
        # Scale to sampling space of cluster image.
        pts2 = pts2 / c_scale
        # Resample as integers for indexing.
        pts2 = np.round(pts2).astype(int)
        # Return nearest cluster values for unassigned points.
        v[v_null] = c_vol[ pts2[v_null, 0], pts2[v_null, 1], pts2[v_null, 2] ]
        # Add in atlas coordinates of unassigned voxels.
        v_pts[v_null, :] = ( pts2[v_null, :] * c_scale ) + c_offset
    # Assemble and return results.
    if map_to_nearest or return_distance or return_voxels:
        v = [v]
    if map_to_nearest:
        v.append(~v_null)
    if return_distance:
        v.append(pdists_min)
    if return_voxels:
        v.append(v_pts)
    return v


def load_deformation_fields(
        df_id='ccf',  # 'ccf' for CCFv3>ARA, or 'ara' for ARA>CCFv3
    ):
    """ Preload deformation fields for mapping between CCFv3 and ARA for repeat usage.
        Note: Deformation fields were computed at 25um resolution.
    """
    assert (df_id == 'ccf') | (df_id == 'ara'), '"id" must be "ccf" or "ara"'
    if df_id == 'ccf':
        df_z = imread(os.path.join(BASE_PATH, r'deformation_fields/CCFv3-to-ARA_25um_deformation_field_z.tif'))
        df_y = imread(os.path.join(BASE_PATH, r'deformation_fields/CCFv3-to-ARA_25um_deformation_field_y.tif'))
        df_x = imread(os.path.join(BASE_PATH, r'deformation_fields/CCFv3-to-ARA_25um_deformation_field_x.tif'))
    elif df_id == 'ara':
        df_z = imread(os.path.join(BASE_PATH, r'deformation_fields/ARA-to-CCFv3_25um_deformation_field_z.tif'))
        df_y = imread(os.path.join(BASE_PATH, r'deformation_fields/ARA-to-CCFv3_25um_deformation_field_y.tif'))
        df_x = imread(os.path.join(BASE_PATH, r'deformation_fields/ARA-to-CCFv3_25um_deformation_field_x.tif'))
    return df_z, df_y, df_x    


def map_coordinates(
        pts,                              # Points to transform, in atlas um resolution [AP, DV, LM].
        df_z=None,                        # Deformation field for z-axis.
        df_y=None,                        # Deformation field for y-axis.
        df_x=None,                        # Deformation field for x-axis.
        df_id='ccf',                      # Which deformation field to load, if not provided.
        scale=np.array([25., 25., 25.]),  # Relative scaling of deformation field (should be 25um).
    ):
    """ Map coordinates using computed deformation fields.  Default is CCFv3>ARA.
        Deformation fields can be pre-loaded and provided for efficiency.
        Note: Deformation fields were computed at 25um resolution.
        If points are outside of deformation field bounds, an error will result.
    """
    # Check for deformation fields, load if necessary.
    if (not df_z) or (not df_y) or (not df_x):
        df_z, df_y, df_x = load_deformation_fields(df_id)
    # Build regular grid interpolator.
    rgi_z = rgi([range(i) for i in df_z.shape], df_z)
    rgi_y = rgi([range(i) for i in df_y.shape], df_y)
    rgi_x = rgi([range(i) for i in df_x.shape], df_x)
    # Scale points to deformation field resolution (25 um).
    pts = pts / scale
    # Map points.
    mapped_pts = np.stack([ rgi_z(pts),
                            rgi_y(pts),
                            rgi_x(pts), ]).T
    # Re-scale points to atlas resolution.
    mapped_pts *= scale
    return mapped_pts


def in_striatum(
        pts,                             # Points to test in atlas um resolution [AP, DV, LM].
        s_vol,                           # Image containing striatum mask.
        s_offset = np.array([0, 0, 0]),  # Offset of s_vol relative to pts, if any (in atlas um resolution).
        s_scale = np.array([1, 1, 1]),   # Scaling of s_vol relative to pts.
    ):
    """ Returns true/false if points are in striatum mask.
    """
    # Init result array.
    tf = np.zeros(pts.shape[0], dtype=bool)
    # Adjust points for offset of striatum mask.
    pts = pts - s_offset
    # Scale points to sampling space of cluster image.
    pts = pts / s_scale
    # Resample as integers for indexing.
    pts = np.round(pts).astype(int)
    # Identify out-of-bounds points (automatically false).
    oob = np.any((pts >= np.array(s_vol.shape)) | (pts < 0), axis=1)
    # Return true/false for in-bounds points only.
    tf[~oob] = s_vol[pts[~oob, 0], pts[~oob, 1], pts[~oob, 2]]
    return tf


def query_striatum_map(
        pts,                       # N x 3 array of points in CCFv3 um space [AP, DV, LM] you want to query.
        cluster_num = 4,           # How many clusters to use (2, 3, 4, or 15), see source figure.
        direction = 'forward',     # "forward" uses cluster maps that are transformed from ARA to CCFv3.
                                   # "reverse" uses original cluster maps in ARA by reverse-transforming query points from CCFv3 to ARA.
                                   # "ara" uses original cluster maps in ARA, and assumes query points are also in ARA um space.
        flip_midline = False,      # Flip any left-side L-M coordinates to right side.
        map_to_nearest = False,    # For pts that are not within defined voxels, return value of nearest voxel.
                                   # This option also returns boolean array indicating which pts are 
                                   # _inside_ voxels (True) vs. outside/mapped (False).
                                   # Note: This option scales poorly with number of points for 'forward', sorry.
        return_distance = False,   # Return euclidean distances from pts to voxel centers.
        return_voxels = False,     # Return coordinates of matched voxels (in atlas um resolution).
        return_instriatum = False, # Return true/false if points are in striatum mask.
                                   # Note: mask is slightly different from cluster map, see readme.
    ):
    """ Primary function to query Hunnicutt et al. (2016) striatum cluster voxels from Fig 5d.
        Note that two different ways to map are possible based on the previously calculated
        transformation between ARA<->CCFv3.  Because this is a non-linear warping, the two
        methods are not perfectly invertible and may result in slight differences for edge cases.
        
        Additional options exist for query points that fall outside of defined cluster voxels.
    """
    # Flip coordinates?
    if flip_midline:
        pts = flip_over_midline(pts)
    # Get parameters based on direction.
    if direction == 'forward':
        # Cluster params.
        c_vol = clusters_ccf[cluster_num]
        c_offset = np.zeros(3)
        c_scale = vx_s['25']
        # Mask params.
        s_vol = mask_ccf
        s_offset = np.zeros(3)
        s_scale = vx_s['25']
    elif (direction == 'reverse') or (direction == 'ara'):
        # Cluster params.
        c_vol = clusters_ara[cluster_num]
        c_offset = vx_t['150']
        c_scale = vx_s['150']
        # Mask params.
        s_vol = mask_ara
        s_offset = np.zeros(3)
        s_scale = vx_s['100']
        if direction == 'reverse':
            # Need to transform points from CCF>ARA.
            pts = map_coordinates(pts, df_id='ccf')
    # Query and return values.
    v = get_cluster_assignments(
            pts,
            c_vol,
            c_offset = c_offset,
            c_scale = c_scale,
            map_to_nearest = map_to_nearest,
            return_distance = return_distance,
            return_voxels = return_voxels,
        )
    # Query striatum mask if requested.
    if return_instriatum:
        tf = in_striatum(
            pts,
            s_vol,
            s_offset = s_offset,
            s_scale = s_scale,
        )
        if map_to_nearest or return_distance or return_voxels:
            v.append(tf)  # Already a list.
        else:
            v = [v, tf]  # Return as list.
    if return_voxels and ((direction == 'reverse') or (direction == 'ara')):
        print('WARNING: With direction=="reverse"/"ara"", returned voxels are in ARA coordinate space!')
    return v


######################## EXAMPLE USAGE BELOW THIS LINE ########################


#### Example data point structure -- N x 3 [A-P, D-V, L-M] um resolution
pts = np.array([
    [ 5842.3, 4182.6, 7914.6, ],
    [ 6123.5, 4823.1, 6752.3, ], 
])

#### Generate random points to test.
pts = np.random.random([500, 3])
pts *= np.array([[7750-3550, 6850-2500, 9550-1850]])
pts += np.array([[3550, 2500, 1850]])
# Map left-side points to right side.
pts = flip_over_midline(pts)

# Visualize points.
viewer.add_points(pts / vx_s['25'], name='random points', scale=vx_s['25'], size=5, face_color='w')

# Get point assignments using 'forward' mapping with 4 cluster dataset.
v_f = query_striatum_map(pts, cluster_num=4, direction='forward')

# Get point assignments using 'reverse' mapping with 4 cluster dataset.
v_r = query_striatum_map(pts, cluster_num=4, direction='reverse')

# How many points returned different values between methods?
print('{:.2f}% of {:d} points differed between the two methods.'.format( sum(v_f != v_r) / len(v_f) * 100, len(v_f)))

# Visualize points showing 'forward' mapped cluster values.
# Dividing and scaling points image resolution to keep napari scrollbar nice.
viewer.layers['cluster[ccf] 04'].visible = True
viewer.add_points(pts / vx_s['25'], name='random points assigned', scale=vx_s['25'], size=5,
    face_color=cluster_colormap['colors'][v_f, :])

# Do forward mapping again, but assign all points based on nearest voxel, even if they are outside cluster bounds.
# This is an extreme example with random points.
v_f2, v_f2_in, v_f2_dist = query_striatum_map(pts, cluster_num=4, direction='forward', map_to_nearest=True, return_distance=True)
viewer.add_points(pts / vx_s['25'], name='random points nearest', scale=vx_s['25'], size=5,
    face_color=cluster_colormap['colors'][v_f2, :])

# Save results to CSV.
csv = pd.DataFrame(pts, columns=['z', 'y', 'x'])
csv['val'] = v_f2
csv['inside'] = v_f2_in
csv['dist'] = v_f2_dist
csv.to_csv(os.path.join(BASE_PATH, r'demo_cluster04.csv'), index=False)


#### Show original cluster voxels as points in ARA space.
# Turn off everything except ARA template.
for item in viewer.layers:
    item.visible = False
ara.visible = True
# Add points.
show_voxels_as_points(clusters_ara[15], c_offset=vx_t['150'], c_scale=vx_s['150'], size=1)


#### Show original cluster voxels as points mapped to CCFv3 space.
# Turn off everything except CCFv3 template.
for item in viewer.layers:
    item.visible = False
ccf.visible = True
# Add points
show_voxels_as_points(clusters_ara[15], c_offset=vx_t['150'], c_scale=vx_s['150'], size=1, tform='ara')


# #### ONE-TIME-USE CODE TO GENERATE STRIATUM CLUSTER VOXELS IN CCFV3 SPACE ####
# 
# from scipy.spatial import Delaunay
# from skimage.io import imsave
# def map_ara_clusters_to_ccf(
#         c_vol,     # Cluster volume to map.
#         c_offset,  # Relative offset of cluster volume.
#         c_scale,   # Scaling of cluster volume (relative to atlas um resolution).
#         ccf,       # Napari image (of CCFv3 25um) to base final output on.
#     ):
#     """ One-time-use routine to create CCFv3 25um-sized volume where voxels are 
#         categorized based on mapping of original ARA-based striatum clusters.
#     """
#     # Set up canvas.
#     canvas = np.zeros(ccf.data.shape).astype(np.uint8)
#     # Quick fix.
#     c_scale = np.array(c_scale)
#     # Get voxel indices in original cluster space.
#     z, y, x = np.where(c_vol)
#     # Get cluster voxel values.
#     v = c_vol[z, y, x]
#     # Arrange voxels indices as N x 3 array and convert to ARA um space.
#     czyx = ( np.array([z, y, x]).T * c_scale ) + c_offset
#     # Initalize "cubes" for each voxel (8 corners are grouped along new axis).
#     cubes = np.tile(czyx, [8, 1, 1]).astype(float)
#     # Apply offsets from voxel centers to get cube vertices.
#     cubes[0, ...] += c_scale * 0.5 * [ 1,  1,  1]
#     cubes[1, ...] += c_scale * 0.5 * [ 1, -1,  1]
#     cubes[2, ...] += c_scale * 0.5 * [ 1, -1, -1]
#     cubes[3, ...] += c_scale * 0.5 * [ 1,  1, -1]
#     cubes[4, ...] += c_scale * 0.5 * [-1,  1,  1]
#     cubes[5, ...] += c_scale * 0.5 * [-1, -1,  1]
#     cubes[6, ...] += c_scale * 0.5 * [-1, -1, -1]
#     cubes[7, ...] += c_scale * 0.5 * [-1,  1, -1]
#     # Reshape cubes to N x 3 list of coordinates.
#     cubes = np.reshape(cubes, [1, -1, 3]).squeeze()
#     # Map cube coordinates from ARA to CCFv3.
#     mapped_cubes = map_coordinates(cubes, df_id='ara')
#     # Scale coordinates to canvas resolution (25 um).
#     mapped_cubes /= ccf.scale
#     # Reshape cube coordinates to regroup individual cube vertices.
#     mapped_cubes = np.reshape(mapped_cubes, [8, -1, 3])
#     # Function to return a small meshgrid that surrounds voxel bounds.
#     def get_bbox_grid(bbox):
#         # Get min/max coordinates of bbox in each dimension.
#         pts_min = np.floor(np.amin(bbox, axis=0)).astype(int)
#         pts_max = np.ceil(np.amax(bbox, axis=0)).astype(int)
#         # Create meshgrid of integer coordinates surrounding bbox in all dimensions.
#         gz, gy, gx = np.meshgrid(*[range(pts_min[i]-1, pts_max[i]+1) for i in range(len(pts_min))])
#         # Reshape meshgrid as N x 3 coordinate list and return.
#         gzyx = np.stack([gz, gy, gx], axis=3).reshape([-1, 3])
#         return gzyx
#     # Loop through each voxel bbox, get 25um-spaced coordinates in vicinity, and test which are "inside" and assign value.
#     # Note: there will be some _literal_ edge cases that appear in two adjacent mapped cubes.
#     # Appears to be ~0.35% of total voxels (3,772 of 1,082,709).  In such cases, second iteration wins, I guess...
#     for i in range(len(v)):
#         bbox = mapped_cubes[:, i, :]
#         gzyx = get_bbox_grid(bbox)
#         idx = Delaunay(bbox).find_simplex(gzyx) >= 0
#         canvas[ gzyx[idx,0], gzyx[idx,1], gzyx[idx,2] ] = v[i]
#     return canvas
# 
# c_res = {}
# for c in clusters_ara:
#     c_res[c] = map_ara_clusters_to_ccf(clusters_ara[c], 
#                                        c_offset=vx_t['150'], c_scale=vx_s['150'], ccf=ccf)
#     imsave(os.path.join(BASE_PATH, r'striatum_clusters_ccfv3/striatum_{:02d}_clusters_ccf.tif'.format(c)), c_res[c], check_contrast=False)
#     viewer.add_image(c_res[c], name='cluster[ccf] {:02d}'.format(c), scale=vx_s['25'],
#                      colormap=cluster_colormap, contrast_limits=[0, 15])                                    			
# 
# # **** TIFs are then opened and resaved in FIJI as OME-TIFFs w/ LZW lossless compression to reduce file size. ****
# 
# ####


# #### ONE-TIME-USE CODE TO GENERATE STRIATUM CLUSTER VOXELS IN CCFV3 SPACE ####
# 
# import numpy as np
# from bg_atlasapi import BrainGlobeAtlas
# from skimage.io import imsave
# 
# # Get CCFv3 at 25um resolution.
# aba = BrainGlobeAtlas('allen_mouse_25um')
# 
# # Generate binary mask for (( dorsal striatum [STRd] + ventral striatum [STRv] ) - olfactory tubercle [OT] )
# strmask = np.logical_xor( np.logical_or( aba.get_structure_mask(aba.structures['STRd']['id']),
#                                          aba.get_structure_mask(aba.structures['STRv']['id']), ),
#                           aba.get_structure_mask(aba.structures['OT']['id'], ))
# 
# # Save mask as TIF.
# imsave(os.path.join(BASE_PATH, r'striatum_masks/striatum_mask_ccf.tif'), strmask.astype('uint8'), check_contrast=False)
# 
# # **** TIFs are then opened and resaved in FIJI as OME-TIFFs w/ LZW lossless compression to reduce file size. ****
# 
# ####