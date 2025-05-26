% Modified from 
%this script shows how everything works together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sunil Bhandari 
%3/17/2017
%Modified Sergio Decenber 2024
%25/12/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [original, movelist, z_slices] = STL_slicer_v3(filename,slice_height,rotate_thetax,rotate_thetay,OUTPUT_read_XLSX)
scale_factor = 1;
triangles = read_stl_FIX_v3(filename,scale_factor);
original = triangles;
triangles = rotate_stl(triangles,'y',90);
num_slices = 100;
[slice_height, z_range] = calculate_slice_height(triangles, num_slices);
[movelist, z_slices] = slice_stl_create_path(triangles, slice_height);

