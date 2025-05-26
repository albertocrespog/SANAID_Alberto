function [Body_Geo,plotData] = Read_data_PLOT_FUS(OUTPUT_read_XLSX,filenameS)

%prefixa = strcat(OUTPUT_read_XLSX.PLOT_flags.fname);
% prefixa = get_fname(OUTPUT_read_XLSX);
% % st0 = strcat('data\');
% % st1 = strcat(st0,prefixa);
% % st1 = strcat(prefixa);
% st1 = filenameS.filename_DATA;
% st2 = strcat('\',OUTPUT_read_XLSX.PLOT_flags.prefix);
% st3 = strcat(st1,st2);

prefixa = get_fname(OUTPUT_read_XLSX);
% st0 = strcat('data\');
% st1 = strcat(st0,prefixa);
% st1 = strcat(prefixa);
st1 = filenameS.filename_DATA;
% st2 = strcat('\',OUTPUT_read_XLSX.PLOT_flags.prefix);
st2 = strcat('\','STORED');
st3 = strcat(st1,st2);
st3A = strcat('\READ_AC.mat');
name = strcat(st3,st3A);

READ_AC = load(name, 'Body_Geo','plotData');

Body_Geo = READ_AC.Body_Geo;
plotData = READ_AC.plotData;
