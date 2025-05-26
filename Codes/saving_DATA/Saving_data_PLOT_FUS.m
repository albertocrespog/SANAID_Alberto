function Saving_data_PLOT_FUS(Body_Geo,plotData,OUTPUT_read_XLSX,filenameS)

%prefixa = strcat(OUTPUT_read_XLSX.PLOT_flags.fname);
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

% Verificar si la carpeta para `name` existe y crearla si no
folder = fileparts(name); % Obtiene solo la ruta de la carpeta
if ~exist(folder, 'dir') % Verifica si la carpeta no existe
    mkdir(folder); % Crea la carpeta
end

save(name, 'Body_Geo','plotData');


