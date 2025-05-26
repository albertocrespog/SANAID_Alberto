function [Storing_GEO_DATA_1,Body_Geo,meshData] = Read_data_Geo(OUTPUT_read_XLSX,filenameS,Geo_tier)

%prefixa = strcat(OUTPUT_read_XLSX.PLOT_flags.fname);
prefixa = get_fname(OUTPUT_read_XLSX);
% st0 = strcat('data\');
% st1 = strcat(st0,prefixa);
% st1 = strcat(prefixa);
st1 = filenameS.filename_DATA;
% st2 = strcat('\',OUTPUT_read_XLSX.PLOT_flags.prefix);
st2 = strcat('\','STORED');
st3 = strcat(st1,st2);
st3A = strcat('\Geo_tier.mat');
name = strcat(st3,st3A);

DATA_GEO = load(name, 'Geo_tier','Body_Geo','meshData');

% Geo_tier = DATA_GEO.Geo_tier;
Body_Geo = DATA_GEO.Body_Geo;
meshData = DATA_GEO.meshData;

Storing_GEO_DATA_1.Geo_tier = Geo_tier;
Storing_GEO_DATA_1.Body_Geo = Body_Geo;
Storing_GEO_DATA_1.meshData = meshData;




