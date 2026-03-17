function [v, varargout] = query_striatum_map(pts, varargin)
% QUERY_STRIATUM_MAP  Query striatum cluster map(s) from Hunnicutt et al. (2016) Fig 5d. 
%   Uses input coordinates in CCFv3 um space [AP, DV, LM].
%   Returns cluster id values (1-15) based on top-to-bottom order of clusters within dendrogram shown in Fig 5c.
%
%   v = QUERY_STRIATUM_MAP(pts) performs default query (4-cluster map, forward mapping).
%
%   Input
%   -------------------------
%   pts: N x 3 array
%       [AP, DV, LM] coordinates to query in CCFv3 um space.
%
%
%   v, varargout = QUERY_STRIATUM_MAP(pts, 'param', 'value') performs query with modified param/value options, see below.
%
%   Optional input parameters
%   -------------------------
%   path: strpath
%       Location of directory containing this script and support files (if not current working directory).
%
%   cluster_num: [15, 4, 3, 2]
%       Which cluster map to use from Fig 5d.
%
%   direction: ['forward', 'reverse', 'ara']
%       'forward' uses cluster maps that are transformed from ARA to CCFv3.
%       'reverse' uses original cluster maps in ARA by reverse-transforming query points from CCFv3 to ARA.
%       'ara' uses original cluster maps in ARA< and assumes query points are also in ARA um space.
%
%   flip_midline: boolean
%       Flip any right-side points over the midline to left side on the M-L axis.
%
%   map_to_nearest: boolean
%       For pts that are not within defined voxels, return value of nearest voxel in varargin.
%
%   return_distance: boolean
%       Return eucldean distances from pts to voxel centers in varargin.
%
%   return_voxels: boolean
%       Return coordinates of matched voxels (in atlas um resolution) in varargin.
%
%   return_instriatum: boolean
%       Return true/false if points are in striatum mask in varargin.
%       Note: mask is slightly different from cluster map, see readme.
%
% v2025.10.13
% muniak@ohsu.edu

% Parse inputs.
p = inputParser;
addRequired(p, 'pts', @isnumeric);
addParameter(p, 'path', pwd, @ischar);
addParameter(p, 'cluster_num', 4, @(x) any(eq(x, [15, 4, 3, 2])));
addParameter(p, 'direction', 'forward', @(x) any(validatestring(x, {'forward', 'reverse', 'ara'})));
addParameter(p, 'flip_midline', false, @islogical);
addParameter(p, 'map_to_nearest', false, @islogical);
addParameter(p, 'return_distance', false, @islogical);
addParameter(p, 'return_voxels', false, @islogical);
addParameter(p, 'return_instriatum', false, @islogical);
parse(p, pts, varargin{:});

% Convenience.
base_path = p.Results.path;
c = p.Results.cluster_num;

% Read in CSVs (though we only need a couple params here).
% Cluster colors (not currently used).
%cluster_colors = csvread(fullfile(base_path, 'csvs', 'cluster_colors.csv'), 1, 1);
% Voxel scales.
tmp = csvread(fullfile(base_path, 'csvs', 'voxel_scales.csv'), 0, 1);
tmp = mat2cell(tmp, [1 3], ones(1, size(tmp, 2)));
tmp(2, :) = cellfun(@transpose, tmp(2, :), 'UniformOutput', false);
vx_s = containers.Map(tmp(1, :), tmp(2, :));
% Voxel translations.
tmp = csvread(fullfile(base_path, 'csvs', 'voxel_translations.csv'), 0, 1);
tmp = mat2cell(tmp, [1 3], ones(1, size(tmp, 2)));
tmp(2, :) = cellfun(@transpose, tmp(2, :), 'UniformOutput', false);
vx_t = containers.Map(tmp(1, :), tmp(2, :));

% Flip coordinates?
if p.Results.flip_midline
    idx = pts(:, 3) < 5700;
    pts(idx, :) = ( pts(idx, :) .* [1 1 -1] ) + [0 0 11400];
end

% Get parameters based on direction.
switch p.Results.direction
    case 'forward'
        % Cluster volume.
        c_vol = read_tif_volume( fullfile(base_path, ...
                                          'striatum_clusters_ccfv3', ...
                                          sprintf('striatum_%02d_clusters_ccf.tif', c)) );
        c_offset = zeros(1, 3);
        c_scale = vx_s(25);
        % Striatum mask.
        if p.Results.return_instriatum
            s_vol = read_tif_volume( fullfile(base_path, ...
                                              'striatum_masks', ...
                                              'striatum_mask_ccf_25um.tif') );
            s_offset = zeros(1, 3);
            s_scale = vx_s(25);
        end
    case {'reverse', 'ara'}
        % Cluster volume.
        c_vol = read_tif_volume( fullfile(base_path, 'striatum_clusters_ara', ...
                                          sprintf('striatum_%02d_clusters_ara.tif', c)) );
        c_offset = vx_t(150);
        c_scale = vx_s(150);
        % Striatum mask.
        if p.Results.return_instriatum
            s_vol = read_tif_volume( fullfile(base_path, ...
                                              'striatum_masks', ...
                                              'striatum_mask_ara-bjh_100um.tif') );
            s_offset = zeros(1, 3);
            s_scale = vx_s(100);
        end
        if strcmp(p.Results.direction, 'reverse')
            % Need to transform points from CCF>ARA.
            pts = map_coordinates(pts);
        end
end

% Test striatum mask first, if requested.
if p.Results.return_instriatum
    % Init result array.
    tf = false(size(pts, 1), 1);
    % Adjust a copy of points for offset of striatum mask.
    pts2 = pts - s_offset;
    % Scale points to sampling space of striatum mask.
    pts2 = pts2  ./ s_scale;
    % Resample as integers for indexing.
    pts2 = round(pts2);
    % python -> matlab...
    pts2 = pts2 + 1;  % +1 to add matlab 1-indexing.
    pts2 = pts2(:, [2 3 1]);  % reorder to matlab image import (YXZ).
    % Identify out-of-bounds points (automatically false).
    oob = any((pts2 > size(s_vol)) | (pts2 <= 0), 2);
    % Return true/false for in-bounds points only.
    tf(~oob) = s_vol(sub2ind(size(s_vol), pts2(~oob, 1), pts2(~oob, 2), pts2 (~oob, 3)));
end

% Init result array.
v = zeros(size(pts, 1), 1, 'uint8');
% Adjust points for offset of cluster image.
pts = pts - c_offset;

% If requested, compute distance between pts and voxel centers before rounding to indices.
if p.Results.return_distance || p.Results.map_to_nearest
    % Get indices of voxel centers.
    [cy, cx, cz] = ind2sub(size(c_vol), find(c_vol > 0));
    % Arrange indices as N x 3 array and scale to atlas resolution
    czyx = ( [cz cy cx] - 1 ) .* c_scale;  % -1 to remove matlab 1-indexing.
    % Find euclidian distances between all points and voxel coordinates.
    pdists = pdist2(pts, czyx, 'euclidean');
    % Find nearest voxel coordinate for each test point, and its euclidian
    % distance.
    [pdists_min, pdists_min_idx] = min(pdists, [], 2);
end

% First, find values for points that are contained in voxels, but not others.
% Scale points to sampling space of cluster image.
pts = pts ./ c_scale;
% Resample as integers for indexing.
pts = round(pts);
% python -> matlab...
pts = pts + 1;  % +1 to add matlab 1-indexing.
pts = pts(:, [2 3 1]);  % reorder to matlab image import (YXZ).
% Identify out-of-bounds points (automatically v==0).
oob = any((pts > size(c_vol)) | (pts <= 0), 2);
% Return cluster values for in-bounds points only.
v(~oob) = c_vol(sub2ind(size(c_vol), pts(~oob, 1), pts(~oob, 2), pts(~oob, 3)));
% matlab -> python...
pts = pts(:, [3 1 2]);  % reorder back to napari indexing (ZYX).
pts = pts - 1;  % -1 to remove matlab 1-indexing.
% Assemble atlas coordinates of assigned voxels.
v_pts = ( pts .* c_scale ) + c_offset;
% Get indices of unassigned points.
v_null = v == 0;
% Set unassigned voxels to nan.
v_pts(v_null, :) = nan;

% If requested, return value of nearest voxel, even if outside of cluster map.
if p.Results.map_to_nearest
    % Get cluster voxel coordinate that was nearest each test point.
    pts2 = czyx(pdists_min_idx, :);
    % Scale to sampling space of cluster image.
    pts2 = pts2 ./ c_scale;
    % Resample as integers for indexing.
    pts2 = round(pts2);
    % Add in atlas coordinates of unassigned voxels.
    v_pts(v_null, :) = ( pts2(v_null, :) .* c_scale ) + c_offset;  % Do this before converting to matlab indexing.
    % python -> matlab...
    pts2 = pts2 + 1;  % +1 to add matlab 1-indexing.
    pts2 = pts2(:, [2 3 1]);  % reorder to matlab image import (YXZ).
    % Return nearest cluster values for unassigned points.
    v(v_null) = c_vol(sub2ind(size(c_vol), pts2(v_null, 1), pts2(v_null, 2), pts2(v_null, 3)));
end

% Assemble optional results.
if p.Results.map_to_nearest || p.Results.return_distance || p.Results.return_voxels || p.Results.return_instriatum
    varargout = {};
end
if p.Results.map_to_nearest
    varargout{end+1} = ~v_null;
end
if p.Results.return_distance
    varargout{end+1} = pdists_min;
end
if p.Results.return_voxels
    varargout{end+1} = v_pts;
end
if p.Results.return_instriatum
    varargout{end+1} = tf;
end

if p.Results.return_voxels && (strcmp(p.Results.direction, 'reverse') || strcmp(p.Results.direction, 'ara'))
    disp('WARNING: With direction=="reverse"/"ara", returned voxels are in ARA coordinate space!');
end


function mapped_pts = map_coordinates(pts)
    % Load previously calculated deformation fields.
    % Note: this matlab function does not contain the versatility of the
    % napari version.
    df_z = read_tif_volume( fullfile(base_path, 'deformation_fields', ...
                                     sprintf('CCFv3-to-ARA_25um_deformation_field_z.tif')) );
    df_y = read_tif_volume( fullfile(base_path, 'deformation_fields', ...
                                     sprintf('CCFv3-to-ARA_25um_deformation_field_y.tif')) );
    df_x = read_tif_volume( fullfile(base_path, 'deformation_fields', ...
                                     sprintf('CCFv3-to-ARA_25um_deformation_field_x.tif')) );
    % Get volume dimensions.
    [sy, sx, sz] = size(df_z);
    % Create mesh grids.
    [gx, gy, gz] = meshgrid(0:sx-1, 0:sy-1, 0:sz-1);  % No need to account for matlab 1-indexing this way.
    % Scale points to deformation field resolution (25 um).
    pts = pts ./ vx_s(25);
    % Map points using grided interpolators.
    mapped_pts = [ interp3(gx, gy, gz, df_z, pts(:, 3), pts(:, 2), pts(:, 1)) ...
                   interp3(gx, gy, gz, df_y, pts(:, 3), pts(:, 2), pts(:, 1)) ...
                   interp3(gx, gy, gz, df_x, pts(:, 3), pts(:, 2), pts(:, 1)) ];
	% Re-scale points to atlas resolution.
    mapped_pts = mapped_pts .* vx_s(25);
end


function tif = read_tif_volume(tif_path)
    if verLessThan('matlab', '9.9')  % tiffreadVolume introduced in R2020b
        im_info = imfinfo(tif_path);
        nz = length(im_info);
        ny = im_info(1).Height;
        nx = im_info(1).Width;
        tif = zeros(ny, nx, nz);
        for z = 1:nz
            tif(:, :, z) = imread(tif_path, 'Index', z);
        end
    else
        tif = tiffreadVolume(tif_path);
    end
end


end