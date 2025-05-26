function Storing_STABILITY_DATA_4D = Saving_Trim_LatD(Trim_ITER_LAT4,Trim_ITER_LAT4B,Trim_ITER_LAT4C,Trim_ITER_LAT4D,conditions_TRIM_lat,OUTPUT_read_XLSX,filenameS)

% Save Mat file
%prefixa = strcat(OUTPUT_read_XLSX.PLOT_flags.fname);
prefixa = get_fname(OUTPUT_read_XLSX);
%st0 = strcat('data\');
%st1 = strcat(st0,prefixa);
%st1 = strcat(prefixa);
st1 = filenameS.filename_DATA;
st2 = strcat('\',OUTPUT_read_XLSX.PLOT_flags.prefix);
st3 = strcat(st1,st2);
st3A = strcat('\Study_Trim_ITER_LAT_D.mat');
name = strcat(st3,st3A);

% Verificar si la carpeta para `name` existe y crearla si no
folder = fileparts(name); % Obtiene solo la ruta de la carpeta
if ~exist(folder, 'dir') % Verifica si la carpeta no existe
    mkdir(folder); % Crea la carpeta
end

save(name, 'Trim_ITER_LAT4','Trim_ITER_LAT4B','Trim_ITER_LAT4C','Trim_ITER_LAT4D','conditions_TRIM_lat')

Storing_STABILITY_DATA_4D.Trim_ITER_LAT4 = Trim_ITER_LAT4;
Storing_STABILITY_DATA_4D.Trim_ITER_LAT4B = Trim_ITER_LAT4B;
Storing_STABILITY_DATA_4D.Trim_ITER_LAT4C = Trim_ITER_LAT4C;
Storing_STABILITY_DATA_4D.Trim_ITER_LAT4D = Trim_ITER_LAT4D;
Storing_STABILITY_DATA_4D.conditions_TRIM_lat = conditions_TRIM_lat;
