%this script shows how everything works together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sunil Bhandari 
%3/17/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear vars; close all;
%triangles = read_binary_stl_file('D638_TypeI.stl');
% triangles = read_ascii_stl('fuselage_A320.stl',1);
triangles = read_ascii_stl('fuselage_EMERGENTIA.stl',1);
original = triangles;
% triangles = orient_stl(triangles,'x');
triangles = rotate_stl(triangles,'y',90);
slice_height = 100;
tic;[movelist, z_slices] = slice_stl_create_path(triangles, slice_height);toc;
%'plotting'
%for i = 1: size(
plot_slices(movelist,z_slices, 0)
%mov = make_movie(movelist, z_slices);
