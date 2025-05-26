% Modified from 
%this script shows how everything works together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sunil Bhandari 
%3/17/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [original, movelist, z_slices] = STL_slicer_v2(filename,slice_height,rotate_thetax,rotate_thetay,OUTPUT_read_XLSX)

% triangles = read_binary_stl_file(fuselage_stl);
% filename1='fuselage_Cessna208_v3_after.stl';

triangles = read_stl_FIX_v2('fuselage_H2_v1_m v1_FUS_ASCII.stl',scale_factor);
%
% if OUTPUT_read_XLSX.Fuselage_flags.STL_reading_test == 1
%     triangles = read_stl_file(filename);
% else OUTPUT_read_XLSX.Fuselage_flags.STL_reading_test == 0
%     triangles = read_ascii_stl(filename,1);
% % elseif OUTPUT_read_XLSX.Fuselage_flags.STL_reading_test == 2
% %     triangles = stlread(filename);
% %     triangles = triangles.ConnectivityList;
% end

original = triangles;

% triangles = orient_stl(triangles,'x');
% Rotates the slices on the y direction 90 degrees so that it creates
% slices in teh y-z plane
triangles = rotate_stl(triangles,'y',rotate_thetay);
triangles = rotate_stl(triangles,'x',rotate_thetax);
% Cuts the slices
[movelist, z_slices] = slice_stl_create_path(triangles, slice_height);
dummy=1;
