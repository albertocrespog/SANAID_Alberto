function Storing_PERFORMANCE_DATA_2 = Saving_data_Performance_var_h_V(Storing_PERFORMANCE_4C,Plot_Options,Restrictions_var_h_V,OUTPUT_read_XLSX,Storing_PROPULSION_DATA,filenameS)
         
Storing_PERFORMANCE_DATA_2.Storing_PERFORMANCE_4C = Storing_PERFORMANCE_4C;
Storing_PERFORMANCE_DATA_2.Plot_Options = Plot_Options;
Storing_PERFORMANCE_DATA_2.Restrictions_var_h_V = Restrictions_var_h_V;
Storing_PERFORMANCE_DATA_2.Storing_PROPULSION_DATA = Storing_PROPULSION_DATA;

% Save Mat file
% prefixa = strcat(OUTPUT_read_XLSX.PLOT_flags.fname);
prefixa = get_fname(OUTPUT_read_XLSX);
% st0 = strcat('data\');
% st1 = strcat(st0,prefixa);
% st1 = strcat(prefixa);
st1 = filenameS.filename_DATA;
st2 = strcat('\',OUTPUT_read_XLSX.PLOT_flags.prefix);
st3 = strcat(st1,st2);
st3A = strcat('\Study_Performance_var_h_var_V.mat');
name = strcat(st3,st3A);

% Verificar si la carpeta para `name` existe y crearla si no
folder = fileparts(name); % Obtiene solo la ruta de la carpeta
if ~exist(folder, 'dir') % Verifica si la carpeta no existe
    mkdir(folder); % Crea la carpeta
end

save(name,'Storing_PERFORMANCE_4C','Plot_Options','Restrictions_var_h_V','Storing_PROPULSION_DATA');