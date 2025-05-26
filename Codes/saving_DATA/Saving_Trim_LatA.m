function Storing_STABILITY_DATA_4A = Saving_Trim_LatA(Trim_ITER_LAT,conditions_TRIM_lat,OUTPUT_read_XLSX,filenameS)

% Save Mat file
%prefixa = strcat(OUTPUT_read_XLSX.PLOT_flags.fname);
prefixa = get_fname(OUTPUT_read_XLSX);
%st0 = strcat('data\');
%st1 = strcat(st0,prefixa);
%st1 = strcat(prefixa);
st1 = filenameS.filename_DATA;
st2 = strcat('\',OUTPUT_read_XLSX.PLOT_flags.prefix);
st3 = strcat(st1,st2);
st3A = strcat('\Study_Trim_ITER_LAT_A.mat');
name = strcat(st3,st3A);

% Verificar si la carpeta para `name` existe y crearla si no
folder = fileparts(name); % Obtiene solo la ruta de la carpeta
if ~exist(folder, 'dir') % Verifica si la carpeta no existe
    mkdir(folder); % Crea la carpeta
end

save(name, 'Trim_ITER_LAT','conditions_TRIM_lat')

Storing_STABILITY_DATA_4A.Trim_ITER_LAT = Trim_ITER_LAT;
Storing_STABILITY_DATA_4A.conditions_TRIM_lat = conditions_TRIM_lat;
