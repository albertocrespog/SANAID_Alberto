function [original, movelist, z_slices_tmp] = Read_data_STL_SLICER(OUTPUT_read_XLSX,filenameS)

%prefixa = strcat(OUTPUT_read_XLSX.PLOT_flags.fname);
prefixa = get_fname(OUTPUT_read_XLSX);
% st0 = strcat('data\');
% st1 = strcat(st0,prefixa);
% st1 = strcat(prefixa);
st1 = filenameS.filename_DATA;
% st2 = strcat('\',OUTPUT_read_XLSX.PLOT_flags.prefix);
st2 = strcat('\','STORED');
st3 = strcat(st1,st2);
st3A = strcat('\STL_SLICER_AC.mat');
name = strcat(st3,st3A);

DATA_STL_SLICER = load(name, 'Geo_tier','Body_Geo','meshData');

original = DATA_STL_SLICER.original;
movelist = DATA_STL_SLICER.movelist;
z_slices_tmp = DATA_STL_SLICER.z_slices_tmp;




