function Storing_STABILITY_DATA_4 = Saving_Trim_Lat(Trim_ITER_LAT,Trim_ITER_LAT2,Trim_ITER_LAT3,conditions_TRIM_lat,OUTPUT_read_XLSX)

% Save Mat file
prefixa = strcat(OUTPUT_read_XLSX.PLOT_flags.fname);
st0 = strcat('data\');
st1 = strcat(st0,prefixa);
st2 = strcat('\',OUTPUT_read_XLSX.PLOT_flags.prefix);
st3 = strcat(st1,st2);
st3A = strcat('\Study_Trim_ITER_LAT.mat');
name = strcat(st3,st3A);
save(name, 'Trim_ITER_LAT','Trim_ITER_LAT2','Trim_ITER_LAT3','conditions_TRIM_lat')

Storing_STABILITY_DATA_4.Trim_ITER_LAT = Trim_ITER_LAT;
Storing_STABILITY_DATA_4.Trim_ITER_LAT2 = Trim_ITER_LAT2;
Storing_STABILITY_DATA_4.Trim_ITER_LAT3 = Trim_ITER_LAT3;
Storing_STABILITY_DATA_4.conditions_TRIM_lat = conditions_TRIM_lat;
